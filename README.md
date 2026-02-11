# RNA-seq Pipeline for Human + Virus Dual-Host Analysis

Snakemake pipeline for paired-end RNA-seq: quality control, alignment to a
combined human + pathogen genome, gene quantification, strandedness-aware
count matrix generation, differential expression (DESeq2), and functional
enrichment (clusterProfiler GO/KEGG).

---

## Prerequisites

### 1. Install Conda (or Mamba)

If you do not have Conda installed, download and install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) or
[Miniforge](https://github.com/conda-forge/miniforge) (which includes Mamba
for faster environment solving):

```bash
# Example: install Miniforge (Linux x86_64)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# Follow the prompts, then restart your shell or run:
source ~/.bashrc
```

After installation, configure strict channel priorities (recommended by
conda-forge and Snakemake):

```bash
conda config --set channel_priority strict
```

### 2. Install Snakemake

Create a dedicated Snakemake environment (recommended) or install into base:

```bash
# Option A: dedicated environment (recommended)
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake

# Option B: install into base
conda install -n base -c conda-forge -c bioconda snakemake
```

### 3. Activate the environment

Before running any `snakemake` command, make sure the environment is active:

```bash
conda activate snakemake   # or 'base' if you installed there
```

---

## Quick Start

```bash
# 1. Clone / copy the project
git clone https://github.com/bruPav/VirusHostRNA.git AdenoGE && cd AdenoGE

# 2. Activate your Snakemake environment
conda activate snakemake

# 3. Place your data (see "Input files" below)

# 4. Edit config.yaml to match your setup

# 5. (Optional) Download reference genomes
snakemake download_references --cores 1

# 6. Run the pipeline
snakemake --use-conda --conda-prefix ~/.snakemake/conda --cores <N>
```

> **Tip:** Using `--conda-prefix ~/.snakemake/conda` stores conda environments
> in a shared location so they are reused across projects and re-runs.

---

## Requirements

| Component   | Minimum | Recommended                    |
|-------------|---------|--------------------------------|
| Snakemake   | >= 7.0  | >= 8.0                         |
| Conda/Mamba | any     | Mamba (faster solving)         |
| CPU cores   | 4       | >= 12                          |
| RAM         | 32 GB   | 48+ GB (STAR indexing)         |
| Disk space  | ~60 GB  | ~100 GB (genome + BAMs + indices) |

All bioinformatics tools (STAR, fastp, samtools, DESeq2, RSeQC,
clusterProfiler, etc.) are installed automatically via conda environments
defined in `envs/`.

---

## Input Files

### 1. FASTQ reads

Place paired-end FASTQ files in `rawdata/` (or the directory set by
`raw_dir` in `config.yaml`). Files must follow this naming convention:

```
rawdata/{sample}_1.fq.gz
rawdata/{sample}_2.fq.gz
```

> **Tip:** If your facility delivers files as `_R1_001.fastq.gz`, rename or
> symlink them, or change `r1_suffix` / `r2_suffix` in `config.yaml`.

Samples are **auto-detected** from filenames -- no need to list them manually.

### 2. Sample metadata (`samples.tsv`)

Create a **tab-separated** file called `samples.tsv` in the project root. See
`samples_template.tsv` for an example.

| Column      | Description |
|-------------|-------------|
| `Sample`    | Must **exactly match** the FASTQ basename (the `{sample}` part). |
| `Condition` | Free-text label for each experimental condition (used for PCA colors, heatmap annotations). |
| `Control`   | `Negative` for control samples, `Experiment` for treated. DESeq2 uses this for the contrast. |

### 3. Reference genomes

The pipeline combines a **human genome** and a **pathogen genome** into a
single reference for alignment. You need four files in the project root:

| File | Default name | Source |
|------|-------------|--------|
| Human genome (gzipped FASTA) | `GRCh38.primary_assembly.genome.fa.gz` | [GENCODE](https://www.gencodegenes.org/human/) |
| Human annotation (gzipped GTF) | `gencode.v49.annotation.gtf.gz` | [GENCODE](https://www.gencodegenes.org/human/) |
| Virus genome (FASTA) | `GCF_006447215.1_ASM644721v1_genomic.fna` | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/) |
| Virus annotation (GTF) | `GCF_006447215.1_ASM644721v1_genomic.gtf` | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/) |

You can auto-download the references (if URLs are set in `config.yaml`):

```bash
snakemake download_references --cores 1
```

For a different pathogen, replace the FASTA/GTF files and update `config.yaml`.

### 4. Combined GTF for strandedness

The rule `make_genes_bed` expects a combined human+pathogen GTF file. By
default `config.yaml` points to `combined_genome.gtf`, which is produced
automatically by the `combine_gtf` step. If you have a hand-curated version
(e.g. with CDS features converted to exon), set `combined_gtf_for_bed` in
`config.yaml` accordingly.

---

## Configuration

All tuneable parameters live in **`config.yaml`**:

```yaml
# Sample sheet
samples_tsv: "samples.tsv"

# Thread limits (capped by --cores at runtime)
threads_star_index: 12
threads_star_align: 8
threads_fastp: 4
threads_fastqc: 2

# STAR read-length parameter (read_length - 1)
sjdb_overhang: 99          # 99 for 100 bp reads, 149 for 150 bp

# Sparse suffix array (set to 2 to halve STAR memory at minor speed cost)
star_genome_sa_sparse_d: 2

# Memory limits (MB)
mem_mb_star_index: 30000
mem_mb_star_align: 8000

# Reference file names, download URLs, FASTQ suffixes, etc.
# ... see config.yaml for the full list
```

---

## Pipeline Overview

```
FASTQ --> fastp (trimming) --> FastQC --> MultiQC
               |
               v
        STAR alignment (human + pathogen combined genome)
               |
               +---> samtools index (BAM indexing)
               |
               +---> RSeQC strandedness inference
               |
               v
        Gene count matrix (strandedness-aware column selection)
               |
               +---> Separate human / viral counts
               |
               +---> Combined     --> DESeq2
               +---> Human counts --> DESeq2 --> Enrichment (GO/KEGG)
               +---> Viral counts --> DESeq2
```

### Rules summary

| Rule | Description |
|------|-------------|
| `download_references` | (Optional) Download human and pathogen genomes from public URLs |
| `unzip_genome` | Decompress human genome FASTA and GTF |
| `prepare_pathogen_genome` | Decompress pathogen files and convert CDS to exon features |
| `combine_genomes` | Concatenate human + pathogen FASTA |
| `combine_gtf` | Concatenate human + pathogen GTF |
| `star_index` | Build STAR genome index on the combined genome |
| `fastp` | Adapter trimming and quality filtering |
| `fastqc` | Per-sample quality reports |
| `multiqc` | Aggregate QC report across all samples |
| `star_align` | Align trimmed reads to the combined genome |
| `samtools_index` | Index sorted BAM files |
| `make_genes_bed` | Convert combined GTF to BED12 for RSeQC |
| `infer_strandedness` | RSeQC `infer_experiment.py` per sample |
| `strandedness_summary` | Summarise strandedness across all samples |
| `count_matrix` | Build gene count matrix (column selected by strandedness) |
| `separate_viral_counts` | Split counts into human and viral matrices |
| `deseq2_analysis` | DESeq2 on combined (human + viral) counts |
| `deseq2_analysis_human` | DESeq2 on human-only counts |
| `deseq2_analysis_viral` | DESeq2 on viral-only counts (sparse-data-safe) |
| `enrichment_analysis` | clusterProfiler ORA + GSEA (GO BP, KEGG) on human DE genes |

### Key output files

| Path | Description |
|------|-------------|
| `results/deg_results_human.tsv` | DESeq2 results for human genes |
| `results/deg_results_significant_human.tsv` | Significant DE genes (padj < 0.05) |
| `results/normalized_counts_human.tsv` | Normalized counts (human) |
| `results/pca_plot_human.png` | PCA of human gene expression |
| `results/volcano_plot_human.png` | Volcano plot (human) |
| `results/heatmap_human.png` | Heatmap of significant genes (human) |
| `results/distance_matrix_human.png` | Sample distance matrix (human) |
| `results/enrichment/` | GO/KEGG ORA and GSEA tables + dot plots |
| `results/deg_results_viral.tsv` | DESeq2 results for pathogen genes |
| `results/deg_results.tsv` | DESeq2 results for combined counts |
| `counts/gene_counts_matrix.tsv` | Combined count matrix |
| `counts/gene_counts_matrix_human.tsv` | Human-only count matrix |
| `counts/gene_counts_matrix_viral.tsv` | Viral-only count matrix |
| `rseqc/strandedness_summary.tsv` | Strandedness summary |
| `multiqc/multiqc_report.html` | Aggregated QC report |

---

## Adapting for a Different Pathogen

1. Replace `pathogen_genome` and `pathogen_gtf` in `config.yaml` with your
   organism's FASTA and GTF files.
2. Update `combined_gtf_for_bed` if you have a custom combined annotation.
3. The `separate_viral_counts` rule identifies pathogen genes by their GTF
   gene IDs and `NC_` chromosome prefixes. If your pathogen uses different
   identifiers, review that rule.

## Adapting for Human-Only (No Pathogen)

If you only have human samples (no co-infection), you can simplify:

1. Set `pathogen_genome` and `pathogen_gtf` to empty or remove the
   `prepare_pathogen_genome`, `combine_genomes`, `combine_gtf`, and
   `separate_viral_counts` rules.
2. Point `star_index` directly to the human genome/GTF.
3. Remove the viral DESeq2 targets from `rule all`.

---

## Conda Environments

| File | Purpose |
|------|---------|
| `envs/ge_analysis.yaml` | Main pipeline: fastp, STAR, samtools, DESeq2, RSeQC, etc. |
| `envs/enrichment.yaml` | clusterProfiler for GO/KEGG enrichment |
| `envs/hipathia.yaml` | HiPathia pathway analysis (optional, currently disabled) |

Environments are created automatically by Snakemake when using `--use-conda`.
The first run will take longer while environments are solved and installed.

---

## Troubleshooting

**STAR runs out of memory during indexing:**
Set `star_genome_sa_sparse_d: 2` in `config.yaml` (halves memory) or reduce
`threads_star_index`. You can also constrain memory with
`--resources mem_mb=32000`.

**biomaRt connection fails:**
DESeq2 maps Ensembl IDs to gene symbols via biomaRt (internet required).
If this fails, results will still contain Ensembl IDs -- the pipeline
does not abort.

**Viral DESeq2 fails ("every gene contains at least one zero"):**
This is expected when viral counts are very sparse. The
`deseq2_analysis_viral.R` script handles this by using `poscounts` size
factors and falling back to log2(normalized + 1) for visualisations when
VST cannot be applied.

**Enrichment analysis finds zero DE genes:**
The enrichment script automatically relaxes `padj` from 0.05 to 0.1 if
fewer than 10 DE genes pass the cutoff.

**Snakemake reports "incomplete files" after an interrupted run:**
Remove the marker files: `rm .snakemake/incomplete/*` and re-run.

---

## Directory Structure

```
AdenoGE/
├── Snakefile                 # Main pipeline definition
├── config.yaml               # All tuneable parameters
├── samples.tsv               # Your sample metadata (create this!)
├── samples_template.tsv      # Template with example entries
├── Samplesheet.tsv           # Sample name mapping (only for standalone HiPathia script, not used by pipeline)
├── envs/                     # Conda environment definitions
│   ├── ge_analysis.yaml
│   ├── enrichment.yaml
│   └── hipathia.yaml
├── scripts/                  # R and Python analysis scripts
│   ├── deseq2_analysis.R
│   ├── deseq2_analysis_viral.R
│   ├── enrichment_analysis_standalone.R
│   ├── gtf2bed.py
│   ├── hipathia_analysis.R
│   └── hipathia_analysis_standalone.R
├── rawdata/                  # Input FASTQ files (you provide)
├── trimmed/                  # fastp output (auto-generated)
├── fastqc/                   # FastQC reports (auto-generated)
├── multiqc/                  # MultiQC report (auto-generated)
├── aligned/                  # STAR BAMs + per-gene counts (auto-generated)
├── star_index/               # STAR genome index (auto-generated)
├── counts/                   # Gene count matrices (auto-generated)
├── results/                  # DESeq2 + enrichment output (auto-generated)
│   └── enrichment/           # GO/KEGG ORA + GSEA tables and plots
├── rseqc/                    # Strandedness reports (auto-generated)
├── refs/                     # BED file for RSeQC (auto-generated)
└── logs/                     # Pipeline logs (auto-generated)
```

---

## License

This pipeline is provided as-is for research purposes.
