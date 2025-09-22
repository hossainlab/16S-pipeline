# 16S rRNA Amplicon Pipeline (QIIME 2 CLI + R/DADA2)

A **production‑grade, fully reproducible** 16S rRNA amplicon pipeline you can run from the terminal on **Ubuntu** or **GitHub Codespaces**.  
Primary track uses **QIIME 2 CLI (Amplicon 2025.7)**; secondary track uses **R/DADA2** to create a parallel set of outputs.  
Everything is automated by `run.sh` and version‑pinned via conda `environment.yml`. QIIME 2 is bootstrapped as its **own, official** environment so it stays conflict‑free.

## Installing Mamba (recommended)

This pipeline requires **Mamba** (a fast drop-in replacement for conda).  
If you don’t already have conda or mamba installed, follow one of the options below.

### Option 1: Install Mambaforge (recommended)

```bash
# Download Mambaforge installer for Linux (x86_64)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

# Run the installer
bash Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge

# Initialize shell (bash example)
$HOME/mambaforge/bin/mamba init bash
source ~/.bashrc
```

```bash
mamba --version
conda --version
```

## Quick start

```bash
# 1) Create the runner env (R/DADA2 + tooling)
mamba env create -f environment.yml
conda activate 16s-pipeline-2025.7

# 2) (One‑time) Bootstrap the official QIIME 2 Amplicon env (2025.7)
bash run.sh --bootstrap-qiime2

# 3) Edit your config
cp config/config.yaml config/myrun.yaml   # adjust primers / trunc lengths / threads
# Place your manifest and metadata in manifests/ and metadata/ (examples provided).

# 4) Run full pipeline (QIIME 2 + DADA2)
bash run.sh -c config/myrun.yaml
```

Outputs land in:
```
outputs/
├── qiime2/   # QIIME 2 artifacts (.qza/.qzv), exports, core metrics
└── dada2/    # Native R/DADA2 ASV table, taxonomy (if refs provided), phyloseq
```

---

## Why two environments?

- **Runner env (`16s-pipeline-2025.7`)** — houses *R/DADA2*, cutadapt, FastQC, MultiQC and CLI utilities for automation. Created from this repo’s `environment.yml`.
- **QIIME 2 env (`qiime2-amplicon-2025.7`)** — installed separately by `run.sh --bootstrap-qiime2` from the **official QIIME 2 Library** YAML for Ubuntu. This preserves QIIME 2’s strict pinning and isolates compiled plugins.

> If you’re in **GitHub Codespaces**, use the provided Dev Container (`.devcontainer/devcontainer.json`): it automatically builds the runner env and then bootstraps QIIME 2.

---

## Repository layout

```
16S-pipeline/
├── environment.yml
├── run.sh
├── README.md
├── .gitignore
├── .github/workflows/ci.yml
├── .devcontainer/devcontainer.json
├── pipeline/
│   └── dada2_track.R
├── config/
│   └── config.yaml
├── manifests/
│   └── manifest.csv       # example (PairedEndFastqManifestPhred33V2)
├── metadata/
│   └── metadata.tsv       # example (QIIME 2 metadata TSV)
├── data/
│   └── raw/               # (optional) local FASTQs
├── ref/
│   ├── classifiers/       # QIIME 2 classifier will auto‑download here
│   └── dada2/             # DADA2 training sets (optional auto‑download; see below)
└── outputs/
    ├── qiime2/
    └── dada2/
```

---

## Configuration (`config/config.yaml`)

Key knobs you’ll likely adjust before the first real run:

```yaml
manifest: "manifests/manifest.csv"          # QIIME 2 manifest (PairedEndFastqManifestPhred33V2)
metadata: "metadata/metadata.tsv"           # sample metadata TSV

# Primers (example V4 515F/806R)
primer_f: "GTGCCAGCMGCCGCGGTAA"
primer_r: "GGACTACHVGGGTWTCTAAT"

# Denoising / trimming
trim_left_f: 0
trim_left_r: 0
trunc_len_f: 240
trunc_len_r: 200

# Diversity
sampling_depth: 10000

# QIIME 2 classifier (will auto‑download + SHA256 verify if missing)
classifier_qza: "ref/classifiers/silva-138-99-nb-classifier.qza"

# Threads
threads: 8

# R/DADA2 options
dada2:
  max_ee_f: 2
  max_ee_r: 2
  trunc_q: 2
  min_len: 50
  pool: "pseudo"           # "pseudo", "true", or "independent"
  tax_train_set: "ref/dada2/silva_nr99_v138.1_train_set.fa.gz"        # optional; see below
  tax_species:    "ref/dada2/silva_species_assignment_v138.1.fa.gz"   # optional
```

> **Tip:** Inspect `outputs/qiime2/demux.qzv` quality plots to pick reasonable `trunc_len_*` and `trim_left_*`. Choose a `sampling_depth` that keeps most samples (use `table.qzv` to assess per‑sample frequency).

---

## Inputs

### QIIME 2 Manifest (Paired‑end, Phred33 V2)
`manifests/manifest.csv`
```csv
sample-id,absolute-filepath,direction
S1,/abs/path/to/S1_R1.fastq.gz,forward
S1,/abs/path/to/S1_R2.fastq.gz,reverse
S2,/abs/path/to/S2_R1.fastq.gz,forward
S2,/abs/path/to/S2_R2.fastq.gz,reverse
```

### Sample metadata (QIIME 2 format)
`metadata/metadata.tsv`
```tsv
#SampleID	Subject	Group
S1	P01	Case
S2	P02	Control
```

---

## What gets produced

**QIIME 2 track** (`outputs/qiime2/`)
- Core artifacts: `demux.qzv`, `table.qzv`, `rep-seqs.qzv`, `taxonomy.qzv`, `taxa-barplot.qzv`
- Phylogeny: `rooted-tree.qza` (plus unrooted)
- Diversity (phylogenetic): `core-metrics/*` (alpha/beta metrics and PCoA)
- Exports (BIOM/TSV/FASTA): `exports/`
- Version manifest: `qiime_info.txt`

**R/DADA2 track** (`outputs/dada2/`)
- `asv_table.tsv`, `seqtab.nochim.rds`, `asv_seqs.fasta`
- `taxonomy.tsv` (if DADA2 training sets provided)
- `phyloseq.rds` (if metadata provided)
- `R_sessionInfo.txt`

---

## Reproducibility & Security

- **Pinned envs**: runner env via `environment.yml`; **QIIME 2** via the official Amplicon 2025.7 YAML (downloaded by URL).
- **Hash checks**: the SILVA classifier (`classifier_qza`) is downloaded and **SHA‑256** verified.
- **Version capture**: `qiime info` and `R sessionInfo()` are saved alongside outputs.
- **Deterministic params**: we fix trimming/denoising options explicitly. DADA2’s algorithm is deterministic for given inputs and parameters.

---

## R/DADA2 reference data

For taxonomy in the native R track you need the DADA2‑formatted **SILVA** training data (e.g., `silva_nr99_v138.1_train_set.fa.gz` and `silva_species_assignment_v138.1.fa.gz`).  
Place them under `ref/dada2/` and point to them in `config.yaml`. If not present, the R pipeline will **skip taxonomy** but still produce ASVs and a `phyloseq` object (without taxonomy).

> You can also train your own region‑specific classifier using QIIME 2’s RESCRIPt plugin and export for R if desired.

---

## GitHub Codespaces (optional)

A ready‑to‑use Dev Container is provided:
```
.devcontainer/devcontainer.json
```
It installs the runner env from `environment.yml` and then bootstraps QIIME 2 automatically.

---

## CI smoke test (optional)

A minimal GitHub Actions workflow `.github/workflows/ci.yml` is included.  
It builds the environment(s) on Ubuntu and runs the pipeline; wire it up to a small test dataset for continuous checks.

---

## Troubleshooting

- **conda/conda‑init**: If `conda run` fails from non‑interactive shells, ensure your shell is conda‑initialized or run via `bash -lc 'conda run -n ENV ...'` inside CI.
- **Manifest paths**: Must be **absolute** for QIIME 2 import. Use `readlink -f` when generating manifests.
- **Sampling depth**: If many samples drop out in diversity, lower `sampling_depth`.
- **Classifier mismatch**: Use a classifier trained for your region (e.g., 515F/806R) for optimal taxonomy; the provided default is full‑length SILVA for broad applicability.
- **Permissions**: After cloning on Linux, mark the runner executable: `chmod +x run.sh`.

---

## License

MIT (see `LICENSE`).

---

## Citation

If you use this pipeline in a publication, please cite:
- QIIME 2, DADA2, SILVA, MAFFT, FastTree, and any reference databases or plugins used.
