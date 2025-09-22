#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dada2)
  library(phyloseq)
  library(ape)
  library(Biostrings)
})

# ---- CLI options ----
option_list <- list(
  make_option(c("--fastq-manifest"), type="character"),
  make_option(c("--metadata"), type="character"),
  make_option(c("--primer-f"), type="character"),
  make_option(c("--primer-r"), type="character"),
  make_option(c("--trunc-len-f"), type="integer", default=0),
  make_option(c("--trunc-len-r"), type="integer", default=0),
  make_option(c("--max-ee-f"), type="double", default=2),
  make_option(c("--max-ee-r"), type="double", default=2),
  make_option(c("--trunc-q"), type="integer", default=2),
  make_option(c("--min-len"), type="integer", default=50),
  make_option(c("--pool"), type="character", default="pseudo"),
  make_option(c("--train-set"), type="character", default=""),
  make_option(c("--species-set"), type="character", default=""),
  make_option(c("--threads"), type="integer", default=8),
  make_option(c("--outdir"), type="character", default="outputs/dada2")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Read QIIME 2 manifest ----
man <- fread(opt$`fastq-manifest`)
# Handle common column name variants
if ("absolute-filepath" %in% names(man)) {
  man[, absolute_file := `absolute-filepath`]
}
if ("filename" %in% names(man) && !("absolute_file" %in% names(man))) {
  man[, absolute_file := normalizePath(filename)]
}
stopifnot(all(c("sample-id","absolute_file","direction") %in% names(man)))

samples <- unique(man$`sample-id`)
fqF <- man[direction=="forward", .(sample=`sample-id`, path=absolute_file)]
fqR <- man[direction=="reverse", .(sample=`sample-id`, path=absolute_file)]
setkey(fqF, sample); setkey(fqR, sample)
stopifnot(identical(fqF$sample, fqR$sample))

# ---- Primer trimming via cutadapt (parity with QIIME 2) ----
message("Trimming primers with cutadapt...")
trimmed_dir <- file.path(opt$outdir, "trimmed")
dir.create(trimmed_dir, showWarnings = FALSE, recursive = TRUE)

run_cutadapt <- function(s, f, r) {
  cmd <- sprintf("cutadapt -j %d -g %s -G %s -o %s/%s_R1.trimmed.fastq.gz -p %s/%s_R2.trimmed.fastq.gz %s %s > %s/%s.cutadapt.log",
                 opt$threads, opt$`primer-f`, opt$`primer-r`,
                 trimmed_dir, s, trimmed_dir, s, f, r, opt$outdir, s)
  system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = FALSE)
}
mapply(run_cutadapt, fqF$sample, fqF$path, fqR$path, SIMPLIFY = FALSE)

# ---- Build paths after trimming ----
filtF <- file.path(trimmed_dir, paste0(fqF$sample, "_R1.trimmed.fastq.gz"))
filtR <- file.path(trimmed_dir, paste0(fqR$sample, "_R2.trimmed.fastq.gz"))
names(filtF) <- fqF$sample; names(filtR) <- fqR$sample

# ---- Filtering ----
filt_dir <- file.path(opt$outdir, "filtered"); dir.create(filt_dir, FALSE, TRUE)
filtFs <- file.path(filt_dir, paste0(names(filtF), "_R1.filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(names(filtR), "_R2.filt.fastq.gz"))

out <- filterAndTrim(filtF, filtFs, filtR, filtRs,
                     truncLen=c(opt$`trunc-len-f`, opt$`trunc-len-r`),
                     maxEE=c(opt$`max-ee-f`, opt$`max-ee-r`),
                     truncQ=opt$`trunc-q`, minLen=opt$`min-len`,
                     rm.phix=TRUE, compress=TRUE, multithread=opt$threads)

# ---- Learn errors & denoise ----
errF <- learnErrors(filtFs, multithread=opt$threads)
errR <- learnErrors(filtRs, multithread=opt$threads)
derepFs <- derepFastq(filtFs); names(derepFs) <- names(filtFs)
derepRs <- derepFastq(filtRs); names(derepRs) <- names(filtRs)
dadaFs <- dada(derepFs, err=errF, multithread=opt$threads, pool=opt$pool)
dadaRs <- dada(derepRs, err=errR, multithread=opt$threads, pool=opt$pool)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=opt$threads)

# ---- Taxonomy (optional) ----
do_tax <- nzchar(opt$`train-set`) && file.exists(opt$`train-set`)
if (do_tax) {
  tax <- assignTaxonomy(colnames(seqtab.nochim), opt$`train-set`, multithread=opt$threads)
  if (nzchar(opt$`species-set`) && file.exists(opt$`species-set`)) {
    tax <- addSpecies(tax, opt$`species-set`)
  }
} else {
  tax <- NULL
  message("No DADA2 training set provided; skipping taxonomy.")
}

# ---- Save outputs ----
saveRDS(seqtab.nochim, file.path(opt$outdir, "seqtab.nochim.rds"))
write.table(t(seqtab.nochim), file.path(opt$outdir, "asv_table.tsv"), sep="\t", quote=FALSE, col.names=NA)

# fasta of ASV sequences
asv_seqs <- Biostrings::DNAStringSet(colnames(seqtab.nochim)); names(asv_seqs) <- paste0("ASV", seq_along(asv_seqs))
Biostrings::writeXStringSet(asv_seqs, filepath=file.path(opt$outdir, "asv_seqs.fasta"))

# taxonomy (if available)
if (!is.null(tax)) {
  write.table(tax, file.path(opt$outdir, "taxonomy.tsv"), sep="\t", quote=FALSE, col.names=NA)
}

# Build a basic phyloseq object when metadata available
if (file.exists(opt$metadata)) {
  meta <- read.table(opt$metadata, header=TRUE, sep="\t", quote="", comment.char="")
  rownames(meta) <- meta[[1]]
  otu <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
  if (!is.null(tax)) {
    tax_tab <- tax_table(as.matrix(tax))
    ps <- phyloseq(otu, tax_tab, sample_data(meta))
  } else {
    ps <- phyloseq(otu, sample_data(meta))
  }
  saveRDS(ps, file.path(opt$outdir, "phyloseq.rds"))
}

# Session info for reproducibility
sink(file.path(opt$outdir, "R_sessionInfo.txt"))
sessionInfo()
sink()
