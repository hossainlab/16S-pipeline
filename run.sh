#!/usr/bin/env bash
set -Eeuo pipefail

# -------- CONFIG / DEFAULTS ---------------------------------------------------
CONFIG=""
Q2_ENV_NAME="qiime2-amplicon-2025.7"
Q2_ENV_YML_URL="https://packages.qiime2.org/amplicon/qiime2-amplicon-2025.7-py310-linux-conda.yml"
# SILVA full-length classifier (QIIME 2 Library) + SHA256
CLASSIFIER_URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
CLASSIFIER_SHA256="c08a1aa4d56b449b511f7215543a43249ae9c54b57491428a7e5548a62613616"

# -------- HELP ----------------------------------------------------------------
usage() {
  cat <<EOF
Usage:
  bash run.sh -c config.yaml            # run both tracks
  bash run.sh --bootstrap-qiime2        # only create QIIME 2 env (idempotent)

Options:
  -c, --config PATH   Path to config YAML (see config/config.yaml)
  -h, --help
EOF
}

# -------- UTIL ----------------------------------------------------------------
timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

log() {
  echo "[$(timestamp)] $*"
}

need() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1"; exit 1; }
}

sha256_verify() {
  local file="$1"; local expect="$2"
  local got
  got=$(sha256sum "$file" | awk '{print $1}')
  if [[ "$got" != "$expect" ]]; then
    echo "ERROR: SHA256 mismatch for $file"
    echo "  got:     $got"
    echo "  expected:$expect"
    exit 1
  fi
}

bootstrap_qiime2() {
  need mamba
  if conda env list | awk '{print $1}' | grep -qx "$Q2_ENV_NAME"; then
    log "QIIME 2 env '$Q2_ENV_NAME' already exists."
    return 0
  fi
  log "Downloading QIIME 2 Amplicon YAML: $Q2_ENV_YML_URL"
  mkdir -p .cache
  curl -L "$Q2_ENV_YML_URL" -o .cache/qiime2-amplicon.yml
  log "Creating env: $Q2_ENV_NAME"
  mamba env create -n "$Q2_ENV_NAME" -f .cache/qiime2-amplicon.yml
  conda run -n "$Q2_ENV_NAME" qiime info | tee outputs/qiime2/qiime_info.txt
}

dl_classifier_if_needed() {
  local path="$1"
  if [[ -f "$path" ]]; then
    log "Classifier exists: $path"
    return 0
  fi
  mkdir -p "$(dirname "$path")"
  log "Downloading SILVA classifier from QIIME 2 Library..."
  curl -L "$CLASSIFIER_URL" -o "$path"
  sha256_verify "$path" "$CLASSIFIER_SHA256"
  log "Classifier downloaded & verified: $path"
}

# -------- ARG PARSE -----------------------------------------------------------
if [[ $# -eq 0 ]]; then usage; exit 1; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--config) CONFIG="$2"; shift 2;;
    --bootstrap-qiime2) bootstrap_qiime2; exit 0;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

# -------- PRECHECKS -----------------------------------------------------------
need mamba; need yq; need jq; need cutadapt; need Rscript
mkdir -p outputs/qiime2 outputs/dada2 logs ref

if [[ ! -f "$CONFIG" ]]; then echo "Config not found: $CONFIG"; exit 1; fi
MANIFEST=$(yq -r '.manifest' "$CONFIG")
METADATA=$(yq -r '.metadata' "$CONFIG")
PRIMER_F=$(yq -r '.primer_f' "$CONFIG")
PRIMER_R=$(yq -r '.primer_r' "$CONFIG")
TRIM_LEFT_F=$(yq -r '.trim_left_f' "$CONFIG")
TRIM_LEFT_R=$(yq -r '.trim_left_r' "$CONFIG")
TRUNC_LEN_F=$(yq -r '.trunc_len_f' "$CONFIG")
TRUNC_LEN_R=$(yq -r '.trunc_len_r' "$CONFIG")
SAMPLING_DEPTH=$(yq -r '.sampling_depth' "$CONFIG")
THREADS=$(yq -r '.threads' "$CONFIG")
CLASSIFIER_QZA=$(yq -r '.classifier_qza' "$CONFIG")
D2_MAXEE_F=$(yq -r '.dada2.max_ee_f' "$CONFIG")
D2_MAXEE_R=$(yq -r '.dada2.max_ee_r' "$CONFIG")
D2_TRUNCQ=$(yq -r '.dada2.trunc_q' "$CONFIG")
D2_MINLEN=$(yq -r '.dada2.min_len' "$CONFIG")
D2_POOL=$(yq -r '.dada2.pool' "$CONFIG")
D2_TRAIN=$(yq -r '.dada2.tax_train_set' "$CONFIG")
D2_SPECIES=$(yq -r '.dada2.tax_species' "$CONFIG")

# Bootstrap QIIME 2 env if missing
bootstrap_qiime2

# Download classifier if needed
dl_classifier_if_needed "$CLASSIFIER_QZA"

# -------- QIIME 2 TRACK -------------------------------------------------------
log "=== QIIME 2: Import demux (manifest: $MANIFEST)"
conda run -n "$Q2_ENV_NAME" qiime tools import   --type "SampleData[PairedEndSequencesWithQuality]"   --input-path "$MANIFEST"   --output-path outputs/qiime2/demux-paired.qza   --input-format PairedEndFastqManifestPhred33V2

conda run -n "$Q2_ENV_NAME" qiime demux summarize   --i-data outputs/qiime2/demux-paired.qza   --o-visualization outputs/qiime2/demux.qzv

log "=== QIIME 2: Primer trimming via q2-cutadapt"
conda run -n "$Q2_ENV_NAME" qiime cutadapt trim-paired   --i-demultiplexed-sequences outputs/qiime2/demux-paired.qza   --p-front-f "$PRIMER_F" --p-front-r "$PRIMER_R"   --p-match-read-wildcards   --p-overlap 5   --o-trimmed-sequences outputs/qiime2/trimmed.qza   --p-cores "$THREADS"

log "=== QIIME 2: DADA2 denoise (paired)"
conda run -n "$Q2_ENV_NAME" qiime dada2 denoise-paired   --i-demultiplexed-seqs outputs/qiime2/trimmed.qza   --p-trunc-len-f "$TRUNC_LEN_F" --p-trunc-len-r "$TRUNC_LEN_R"   --p-trim-left-f "$TRIM_LEFT_F" --p-trim-left-r "$TRIM_LEFT_R"   --p-n-threads "$THREADS"   --o-table outputs/qiime2/table.qza   --o-representative-sequences outputs/qiime2/rep-seqs.qza   --o-denoising-stats outputs/qiime2/denoising-stats.qza

conda run -n "$Q2_ENV_NAME" qiime feature-table summarize   --i-table outputs/qiime2/table.qza   --m-sample-metadata-file "$METADATA"   --o-visualization outputs/qiime2/table.qzv

conda run -n "$Q2_ENV_NAME" qiime feature-table tabulate-seqs   --i-data outputs/qiime2/rep-seqs.qza   --o-visualization outputs/qiime2/rep-seqs.qzv

log "=== QIIME 2: Phylogeny"
conda run -n "$Q2_ENV_NAME" qiime phylogeny align-to-tree-mafft-fasttree   --i-sequences outputs/qiime2/rep-seqs.qza   --o-alignment outputs/qiime2/aligned-rep-seqs.qza   --o-masked-alignment outputs/qiime2/masked-aligned-rep-seqs.qza   --o-tree outputs/qiime2/unrooted-tree.qza   --o-rooted-tree outputs/qiime2/rooted-tree.qza

log "=== QIIME 2: Taxonomy (sklearn) with SILVA 138 full-length classifier"
conda run -n "$Q2_ENV_NAME" qiime feature-classifier classify-sklearn   --i-classifier "$CLASSIFIER_QZA"   --i-reads outputs/qiime2/rep-seqs.qza   --o-classification outputs/qiime2/taxonomy.qza

conda run -n "$Q2_ENV_NAME" qiime metadata tabulate   --m-input-file outputs/qiime2/taxonomy.qza   --o-visualization outputs/qiime2/taxonomy.qzv

conda run -n "$Q2_ENV_NAME" qiime taxa barplot   --i-table outputs/qiime2/table.qza   --i-taxonomy outputs/qiime2/taxonomy.qza   --m-metadata-file "$METADATA"   --o-visualization outputs/qiime2/taxa-barplot.qzv

log "=== QIIME 2: Core metrics (phylogenetic)"
conda run -n "$Q2_ENV_NAME" qiime diversity core-metrics-phylogenetic   --i-phylogeny outputs/qiime2/rooted-tree.qza   --i-table outputs/qiime2/table.qza   --p-sampling-depth "$SAMPLING_DEPTH"   --m-metadata-file "$METADATA"   --p-n-jobs-or-threads "$THREADS"   --output-dir outputs/qiime2/core-metrics

log "=== QIIME 2: Exports (BIOM/TSV/FASTA)"
mkdir -p outputs/qiime2/exports
conda run -n "$Q2_ENV_NAME" qiime tools export --input-path outputs/qiime2/table.qza --output-path outputs/qiime2/exports
conda run -n "$Q2_ENV_NAME" qiime tools export --input-path outputs/qiime2/taxonomy.qza --output-path outputs/qiime2/exports
conda run -n "$Q2_ENV_NAME" qiime tools export --input-path outputs/qiime2/rep-seqs.qza --output-path outputs/qiime2/exports

# Capture versions for reproducibility
conda run -n "$Q2_ENV_NAME" qiime info > outputs/qiime2/qiime_info.txt

# -------- R/DADA2 TRACK -------------------------------------------------------
log "=== R/DADA2: running native pipeline"
# If DADA2 training sets are missing, taxonomy will be skipped by the script.
Rscript pipeline/dada2_track.R   --fastq-manifest "$MANIFEST"   --metadata "$METADATA"   --primer-f "$PRIMER_F"   --primer-r "$PRIMER_R"   --trunc-len-f "$TRUNC_LEN_F"   --trunc-len-r "$TRUNC_LEN_R"   --max-ee-f "$D2_MAXEE_F"   --max-ee-r "$D2_MAXEE_R"   --trunc-q "$D2_TRUNCQ"   --min-len "$D2_MINLEN"   --pool "$D2_POOL"   --train-set "$D2_TRAIN"   --species-set "$D2_SPECIES"   --threads "$THREADS"   --outdir outputs/dada2

log "DONE. See outputs/qiime2 and outputs/dada2"
