# AdenoGE — RNA-seq Pipeline for Human + Adenovirus Dual-Host Analysis

Snakemake pipeline for paired-end RNA-seq: quality control, alignment to a
combined human + pathogen genome, gene quantification, differential expression
(DESeq2), and functional enrichment (clusterProfiler).

---

## Quick Start

```bash
# 1. Clone / copy the project
git clone <repo-url> AdenoGE && cd AdenoGE

# 2. Install Snakemake (if not already available)
conda install -n base -c conda-forge -c bioconda snakemake

# 3. Place your data (see "Input files" below)

# 4. Edit config.yaml to match your setup

# 5. Run the pipeline
snakemake --use-conda --cores <N>
```

---

## Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| Snakemake | >= 7.0 | >= 8.0 |
| Conda / Mamba | any | Mamba (faster solving) |
| CPU cores | 4 | >= 12 |
| RAM | 32 GB | 48+ GB (STAR indexing) |
| Disk space | ~60 GB | ~100 GB (genome + BAMs + indices) |

All bioinformatics tools (STAR, fastp, DESeq2, etc.) are installed
automatically via conda environments defined in `envs/`.

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

Samples are **auto-detected** from filenames — no need to list them manually.

### 2. Sample metadata (`samples.tsv`)

Create a **tab-separated** file called `samples.tsv` in the project root. See
`samples_template.tsv` for an example.

| Column | Description |
|--------|-------------|
| `Sample` | Must **exactly match** the FASTQ basename (the `{sample}` part). |
| `Condition` | Free-text label for each experimental condition (used for PCA colors, heatmap annotations). |
| `Control` | `Negative` for control samples, `Experiment` for treated. DESeq2 uses this for the contrast. |
| `Exp_group_1` | Integer grouping variable (e.g., cell line). |
| `Exp_group_2` | Integer grouping variable (e.g., time point). |

### 3. Reference genomes

The pipeline combines a **human genome** and a **pathogen genome** into a
single reference for alignment. You need four files in the project root:

| File | Default name | Source |
|------|-------------|--------|
| Human genome (gzipped FASTA) | `GRCh38.primary_assembly.genome.fa.gz` | [GENCODE](https://www.gencodegenes.org/human/) |
| Human annotation (gzipped GTF) | `gencode.v49.annotation.gtf.gz` | [GENCODE](https://www.gencodegenes.org/human/) |
| Pathogen genome (FASTA) | `GCF_006447215.1_ASM644721v1_genomic.fna` | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/) |
| Pathogen annotation (GTF) | `GCF_006447215.1_ASM644721v1_genomic.gtf` | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/) |

You can also auto-download the human references:

```bash
snakemake download_references --cores 1
```

For a different pathogen, replace the FASTA/GTF files and update `config.yaml`.

### 4. Combined GTF for strandedness

The rule `make_genes_bed` expects a combined human+pathogen GTF file
(`combined_genome.CDSasExon.gtf` by default). If you do not have this file,
the pipeline will create `combined_genome.gtf` during the `combine_gtf` step;
point `combined_gtf_for_bed` in `config.yaml` to `combined_genome.gtf` instead.

---

## Configuration

All tuneable parameters live in **`config.yaml`**:

```yaml
# Sample sheet path
samples_tsv: "samples.tsv"

# Thread limits (capped by --cores at runtime)
threads_star_index: 12
threads_star_align: 8
threads_fastp: 4
threads_fastqc: 2

# STAR read-length parameter (read_length - 1)
sjdb_overhang: 99          # 99 for 100 bp reads, 149 for 150 bp

# Reference file names
human_genome_gz: "GRCh38.primary_assembly.genome.fa.gz"
pathogen_genome: "your_pathogen.fna"
# ... see config.yaml for the full list
```

---

## Pipeline Overview

```
FASTQ ──► fastp (trimming) ──► FastQC ──► MultiQC
                │
                ▼
         STAR alignment (human + pathogen combined genome)
                │
                ├──► RSeQC strandedness inference
                │
                ▼
         Gene count matrix (strandedness-aware column selection)
                │
                ├──► Human counts ──► DESeq2 ──► Enrichment (GO/KEGG)
                ├──► Viral counts ──► DESeq2
                └──► Combined     ──► DESeq2
```

### Key output files

| Path | Description |
|------|-------------|
| `results/deg_results_human.tsv` | DESeq2 results for human genes |
| `results/deg_results_significant_human.tsv` | Significant DE genes (padj < 0.05) |
| `results/pca_plot_human.png` | PCA of human gene expression |
| `results/volcano_plot_human.png` | Volcano plot |
| `results/heatmap_human.png` | Heatmap of significant genes |
| `results/enrichment/` | GO/KEGG ORA and GSEA tables + plots |
| `results/deg_results_viral.tsv` | DESeq2 results for pathogen genes |
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
   `prepare_adenovirus_genome`, `combine_genomes`, `combine_gtf`, and
   `separate_viral_counts` rules.
2. Point `star_index` directly to the human genome/GTF.
3. Remove the viral DESeq2 targets from `rule all`.

---

## Conda Environments

| File | Purpose |
|------|---------|
| `envs/ge_analysis.yaml` | Main pipeline: fastp, STAR, samtools, DESeq2, etc. |
| `envs/enrichment.yaml` | clusterProfiler for GO/KEGG enrichment |
| `envs/hipathia.yaml` | HiPathia pathway analysis (optional) |

Environments are created automatically by Snakemake when using `--use-conda`.

---

## Troubleshooting

**STAR runs out of memory during indexing:**
Reduce `threads_star_index` in `config.yaml` or run with
`--resources mem_mb=32000` to constrain Snakemake.

**biomaRt connection fails:**
DESeq2 maps Ensembl IDs to gene symbols via biomaRt (internet required).
If this fails, results will still contain Ensembl IDs — the pipeline
does not abort.

**Enrichment analysis finds zero DE genes:**
The enrichment script automatically relaxes `padj` from 0.05 to 0.1 if
fewer than 10 DE genes pass the cutoff.

---

## Directory Structure

```
AdenoGE/
├── Snakefile                 # Main pipeline
├── config.yaml               # All tuneable parameters
├── samples.tsv               # Your sample metadata (create this!)
├── samples_template.tsv      # Template with example entries
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
├── aligned/                  # STAR BAMs + counts (auto-generated)
├── star_index/               # STAR index (auto-generated)
├── counts/                   # Gene count matrices (auto-generated)
├── results/                  # DESeq2 + enrichment output (auto-generated)
├── rseqc/                    # Strandedness reports (auto-generated)
├── refs/                     # BED file for RSeQC (auto-generated)
└── logs/                     # Pipeline logs (auto-generated)
```

---

## License

This pipeline is provided as-is for research purposes.
