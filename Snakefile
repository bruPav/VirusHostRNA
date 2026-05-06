import glob
import os

# =============================================================================
# Load configuration
# =============================================================================
configfile: "config.yaml"

# --- Directories (from config, with sensible defaults) -----------------------
RAW_DIR     = config.get("raw_dir", "rawdata")
TRIM_DIR    = config.get("trim_dir", "trimmed")
FASTQC_DIR  = config.get("fastqc_dir", "fastqc")
MULTIQC_DIR = config.get("multiqc_dir", "multiqc")
ALIGN_DIR   = config.get("align_dir", "aligned")
INDEX_DIR   = config.get("index_dir", "star_index")
RESULTS_DIR = config.get("results_dir", "results")
RSEQC_DIR   = config.get("rseqc_dir", "rseqc")
COUNTS_DIR  = config.get("counts_dir", "counts")
REFS_DIR    = config.get("refs_dir", "refs")
LOGS_DIR    = config.get("logs_dir", "logs")

# --- Reference files ----------------------------------------------------------
HUMAN_GENOME_GZ = config.get("human_genome_gz", "GRCh38.primary_assembly.genome.fa.gz")
HUMAN_GTF_GZ    = config.get("human_gtf_gz", "gencode.v49.annotation.gtf.gz")
PATHOGEN_GENOME = config.get("pathogen_genome", "GCF_006447215.1_ASM644721v1_genomic.fna")
PATHOGEN_GTF    = config.get("pathogen_gtf", "GCF_006447215.1_ASM644721v1_genomic.gtf")
COMBINED_GTF_BED = config.get("combined_gtf_for_bed", "combined_genome.gtf")

# Derived names (uncompressed)
HUMAN_GENOME = HUMAN_GENOME_GZ.replace(".gz", "")
HUMAN_GTF    = HUMAN_GTF_GZ.replace(".gz", "")

# Pathogen FA: strip any known extension (.fa, .fna, .fna.gz, .fa.gz, etc.)
# and append _prepared.fa to avoid a circular dependency (same path for input/output).
_pathogen_basename = os.path.basename(PATHOGEN_GENOME)
for ext in ['.fna.gz', '.fa.gz', '.fna', '.fa']:
    if _pathogen_basename.endswith(ext):
        _pathogen_basename = _pathogen_basename[:-len(ext)]
        break
PATHOGEN_FA = _pathogen_basename + "_prepared.fa"
PATHOGEN_GTF_PREPARED = os.path.splitext(PATHOGEN_GTF)[0] + "_prepared.gtf"

COMBINED_GENOME = "combined_genome.fa"
COMBINED_GTF    = "combined_genome.gtf"
GENES_BED       = os.path.join(REFS_DIR, "genes.bed")

# --- FASTQ naming convention -------------------------------------------------
R1_SUFFIX = config.get("r1_suffix", "_1.fq.gz")
R2_SUFFIX = config.get("r2_suffix", "_2.fq.gz")

LIBRARY_TYPE = config.get("library_type", "PE")
IS_PAIRED    = LIBRARY_TYPE == "PE"
SE_SUFFIX    = config.get("se_suffix", ".fq.gz")

# --- Samples (auto-detected from raw_dir) ------------------------------------
if IS_PAIRED:
    SAMPLES = sorted(
        set(
            os.path.basename(f).replace(R1_SUFFIX, "")
            for f in glob.glob(os.path.join(RAW_DIR, f"*{R1_SUFFIX}"))
        )
    )
else:
    SAMPLES = sorted(
        set(
            os.path.basename(f).replace(SE_SUFFIX, "")
            for f in glob.glob(os.path.join(RAW_DIR, f"*{SE_SUFFIX}"))
        )
    )

if not SAMPLES:
    suffix = R1_SUFFIX if IS_PAIRED else SE_SUFFIX
    raise SystemExit(
        f"ERROR: No FASTQ files found in '{RAW_DIR}' matching '*{suffix}'.\n"
        f"Place your {'paired-end' if IS_PAIRED else 'single-end'} reads there and re-run."
    )

# --- Sample metadata ----------------------------------------------------------
SAMPLES_TSV = config.get("samples_tsv", "samples.tsv")

# --- STAR parameters ----------------------------------------------------------
SJDB_OVERHANG = config.get("sjdb_overhang", 99)

# --- Thread limits (capped by workflow.cores at runtime) ----------------------
THREADS_STAR_INDEX = config.get("threads_star_index", 12)
THREADS_STAR_ALIGN = config.get("threads_star_align", 8)
THREADS_FASTP      = config.get("threads_fastp", 4)
THREADS_FASTQC     = config.get("threads_fastqc", 2)

# --- Memory limits (MB) ------------------------------------------------------
MEM_MB_STAR_INDEX = config.get("mem_mb_star_index", 30000)
MEM_MB_STAR_ALIGN = config.get("mem_mb_star_align", 8000)

# --- STAR sparse SA (reduces memory for genome generation) -------------------
STAR_SA_SPARSE_D = config.get("star_genome_sa_sparse_d", 2)

# =============================================================================
# Pre-flight validation — catch misconfigurations before any rule runs
# =============================================================================
def _validate_config():
    """Check config, input files, and metadata for common errors."""
    errors = []

    # --- Library type ---
    lib_type = config.get("library_type", "PE")
    if lib_type not in ("PE", "SE"):
        errors.append(f"library_type must be 'PE' or 'SE', got '{lib_type}'.")
    if lib_type == "SE" and not config.get("se_suffix"):
        errors.append("se_suffix must be set when library_type is 'SE'.")

    # --- Integer configs must be positive ---
    for key in ["threads_star_index", "threads_star_align", "threads_fastp",
                "threads_fastqc", "mem_mb_star_index", "mem_mb_star_align",
                "sjdb_overhang"]:
        v = config.get(key)
        if v is not None and (not isinstance(v, (int, float)) or v <= 0):
            errors.append(f"{key} must be a positive number, got {v}.")

    # --- pathogen_gene_prefixes must be a list ---
    pfx = config.get("pathogen_gene_prefixes")
    if pfx is not None and not isinstance(pfx, list):
        errors.append(f"pathogen_gene_prefixes must be a list, got {type(pfx).__name__}.")

    # --- Reference files must exist ---
    for key, path in [("human_genome_gz", HUMAN_GENOME_GZ),
                      ("human_gtf_gz",    HUMAN_GTF_GZ),
                      ("pathogen_genome",  PATHOGEN_GENOME),
                      ("pathogen_gtf",     PATHOGEN_GTF)]:
        if not os.path.isfile(path):
            errors.append(
                f"Reference file '{path}' ({key}) not found.\n"
                f"  Place it in the project root or run: snakemake download_references --cores 1"
            )

    # --- Sample metadata ---
    if not os.path.isfile(SAMPLES_TSV):
        errors.append(
            f"Sample metadata file '{SAMPLES_TSV}' not found.\n"
            f"  Copy samples_template.tsv to {SAMPLES_TSV} and edit it."
        )
    else:
        with open(SAMPLES_TSV) as f:
            header = f.readline().strip().split("\t")
        required_cols = ["Sample"]
        missing = [c for c in required_cols if c not in header]
        if missing:
            errors.append(
                f"Missing required column(s) {missing} in '{SAMPLES_TSV}'.\n"
                f"  Expected columns: {required_cols}"
            )
        if "Sample" in header:
            metadata_samples = set()
            with open(SAMPLES_TSV) as f:
                next(f)
                for line in f:
                    parts = line.strip().split("\t")
                    if parts:
                        metadata_samples.add(parts[0])
            missing_from_meta = set(SAMPLES) - metadata_samples
            missing_from_data = metadata_samples - set(SAMPLES)
            if missing_from_meta:
                errors.append(
                    f"Samples in FASTQ but not in '{SAMPLES_TSV}': {sorted(missing_from_meta)}"
                )
            if missing_from_data:
                errors.append(
                    f"Samples in '{SAMPLES_TSV}' but no FASTQ found: {sorted(missing_from_data)}"
                )

    # --- Output directories should be writable ---
    for d in [TRIM_DIR, FASTQC_DIR, MULTIQC_DIR, ALIGN_DIR, RESULTS_DIR,
              RSEQC_DIR, COUNTS_DIR, REFS_DIR, LOGS_DIR]:
        try:
            os.makedirs(d, exist_ok=True)
        except OSError as e:
            errors.append(f"Cannot create/write to '{d}': {e}")

    if errors:
        msg = "\n  - ".join(errors)
        raise SystemExit(
            f"Pre-flight validation failed with {len(errors)} error(s):\n  - {msg}"
        )

_validate_config()

# =============================================================================
# Target rule
# =============================================================================
def _rule_all_inputs():
    """Build the full list of target outputs, conditional on library type."""
    ins = []

    # Trimmed reads
    if IS_PAIRED:
        ins.extend(expand(f"{TRIM_DIR}/{{sample}}_1.trim.fq.gz", sample=SAMPLES))
    else:
        ins.extend(expand(f"{TRIM_DIR}/{{sample}}.trim.fq.gz", sample=SAMPLES))

    # FastQC
    if IS_PAIRED:
        ins.extend(expand(f"{FASTQC_DIR}/{{sample}}_1.trim_fastqc.html", sample=SAMPLES))
    else:
        ins.extend(expand(f"{FASTQC_DIR}/{{sample}}.trim_fastqc.html", sample=SAMPLES))

    # Sorted BAMs, BAM indices, and gene counts
    ins.extend(expand(f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES))
    ins.extend(expand(f"{ALIGN_DIR}/{{sample}}_BAMIndex.bai", sample=SAMPLES))
    ins.extend(expand(f"{ALIGN_DIR}/{{sample}}_ReadsPerGene.out.tab", sample=SAMPLES))
    ins.append(f"{MULTIQC_DIR}/multiqc_report.html")
    # Strandedness reports
    ins.extend(expand(f"{RSEQC_DIR}/{{sample}}_infer_experiment.txt", sample=SAMPLES))
    ins.append(f"{RSEQC_DIR}/strandedness_summary.tsv")
    ins.append(f"{COUNTS_DIR}/gene_counts_matrix.tsv")
    ins.append(f"{COUNTS_DIR}/gene_counts_matrix_human.tsv")
    ins.append(f"{COUNTS_DIR}/gene_counts_matrix_viral.tsv")
    # Combined DESeq2
    ins.append(f"{RESULTS_DIR}/deg_results.tsv")
    ins.append(f"{RESULTS_DIR}/deg_results_significant.tsv")
    ins.append(f"{RESULTS_DIR}/normalized_counts.tsv")
    ins.append(f"{RESULTS_DIR}/pca_plot.png")
    ins.append(f"{RESULTS_DIR}/distance_matrix.png")
    ins.append(f"{RESULTS_DIR}/volcano_plot.png")
    ins.append(f"{RESULTS_DIR}/heatmap.png")
    # Human DESeq2
    ins.append(f"{RESULTS_DIR}/deg_results_human.tsv")
    ins.append(f"{RESULTS_DIR}/deg_results_significant_human.tsv")
    ins.append(f"{RESULTS_DIR}/normalized_counts_human.tsv")
    ins.append(f"{RESULTS_DIR}/pca_plot_human.png")
    ins.append(f"{RESULTS_DIR}/distance_matrix_human.png")
    ins.append(f"{RESULTS_DIR}/volcano_plot_human.png")
    ins.append(f"{RESULTS_DIR}/heatmap_human.png")
    # Viral DESeq2
    ins.append(f"{RESULTS_DIR}/deg_results_viral.tsv")
    ins.append(f"{RESULTS_DIR}/deg_results_significant_viral.tsv")
    ins.append(f"{RESULTS_DIR}/normalized_counts_viral.tsv")
    ins.append(f"{RESULTS_DIR}/pca_plot_viral.png")
    ins.append(f"{RESULTS_DIR}/distance_matrix_viral.png")
    ins.append(f"{RESULTS_DIR}/volcano_plot_viral.png")
    ins.append(f"{RESULTS_DIR}/heatmap_viral.png")
    # Enrichment
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_KEGG_GSEA.tsv")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_GOBP_GSEA.tsv")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_KEGG_GSEA_dotplot.png")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_GOBP_GSEA_dotplot.png")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_KEGG_ORA.tsv")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_GOBP_ORA.tsv")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_KEGG_ORA_dotplot.png")
    ins.append(f"{RESULTS_DIR}/enrichment/enrich_GOBP_ORA_dotplot.png")
    # Summary report
    ins.append(f"{RESULTS_DIR}/report.html")
    # Version log
    ins.append(f"{RESULTS_DIR}/versions.txt")
    return ins

rule all:
    input:
        _rule_all_inputs()

# =============================================================================
# Download reference genomes (optional convenience rule)
# =============================================================================
def _ref_download_targets():
    """Only declare outputs for references that have a download URL configured."""
    targets = []
    if config.get("human_genome_url"):
        targets.append(HUMAN_GENOME_GZ)
    if config.get("human_gtf_url"):
        targets.append(HUMAN_GTF_GZ)
    if config.get("pathogen_genome_url"):
        targets.append(PATHOGEN_GENOME)
    if config.get("pathogen_gtf_url"):
        targets.append(PATHOGEN_GTF)
    return targets

rule download_references:
    """Download reference genomes from public repositories.
    Run explicitly: snakemake download_references --cores 1
    Only files whose URL is set (non-empty) in config.yaml are downloaded.
    """
    output:
        _ref_download_targets()
    params:
        url_map = [
            (HUMAN_GENOME_GZ,  config.get("human_genome_url", "")),
            (HUMAN_GTF_GZ,     config.get("human_gtf_url", "")),
            (PATHOGEN_GENOME,   config.get("pathogen_genome_url", "")),
            (PATHOGEN_GTF,      config.get("pathogen_gtf_url", "")),
        ]
    run:
        import subprocess
        for target, url in params.url_map:
            if url:
                print(f"Downloading {target} from {url} ...")
                subprocess.check_call(["wget", "-c", "-O", target, url])
            else:
                print(f"Skipping {target} (no URL configured in config.yaml)")

# =============================================================================
# Genome preparation
# =============================================================================
rule unzip_genome:
    input:
        fa_gz  = HUMAN_GENOME_GZ,
        gtf_gz = HUMAN_GTF_GZ
    output:
        fa  = HUMAN_GENOME,
        gtf = HUMAN_GTF
    shell:
        """
        gunzip -c {input.fa_gz} > {output.fa}
        gunzip -c {input.gtf_gz} > {output.gtf}
        """

rule prepare_pathogen_genome:
    """Decompress pathogen FASTA/GTF if gzipped and convert CDS→exon.

    Many pathogen GTFs (e.g. NCBI RefSeq) only annotate CDS features and
    lack exon features.  STAR --quantMode GeneCounts only counts reads that
    overlap "exon" features, so without the CDS→exon conversion the
    pathogen genes would get zero counts.
    """
    input:
        fa  = PATHOGEN_GENOME,
        gtf = PATHOGEN_GTF
    output:
        fa_prepared  = PATHOGEN_FA,
        gtf_prepared = PATHOGEN_GTF_PREPARED
    shell:
        """
        # Decompress FASTA if gzipped, otherwise copy
        if file {input.fa} | grep -q "gzip"; then
            gunzip -c {input.fa} > {output.fa_prepared}
        else
            cp {input.fa} {output.fa_prepared}
        fi

        # Decompress GTF if gzipped, then convert CDS features to exon.
        if file {input.gtf} | grep -q "gzip"; then
            gunzip -c {input.gtf} | awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{if($3=="CDS") $3="exon"; print}}' > {output.gtf_prepared}
        else
            awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{if($3=="CDS") $3="exon"; print}}' {input.gtf} > {output.gtf_prepared}
        fi
        """

rule combine_genomes:
    input:
        human_fa     = HUMAN_GENOME,
        pathogen_fa  = PATHOGEN_FA
    output:
        combined = COMBINED_GENOME
    shell:
        """
        cat {input.human_fa} {input.pathogen_fa} > {output.combined}
        """

rule combine_gtf:
    input:
        human_gtf    = HUMAN_GTF,
        pathogen_gtf = PATHOGEN_GTF_PREPARED
    output:
        combined = COMBINED_GTF
    shell:
        """
        cat {input.human_gtf} {input.pathogen_gtf} > {output.combined}
        """

# =============================================================================
# STAR genome index
# =============================================================================
rule star_index:
    input:
        fa  = COMBINED_GENOME,
        gtf = COMBINED_GTF
    output:
        directory(INDEX_DIR)
    threads: THREADS_STAR_INDEX
    resources:
        mem_mb = MEM_MB_STAR_INDEX
    conda: "envs/ge_analysis.yaml"
    params:
        sa_sparse = STAR_SA_SPARSE_D
    shell:
        """
        mkdir -p {output}
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fa} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {SJDB_OVERHANG} \
             --genomeSAsparseD {params.sa_sparse}
        """

# =============================================================================
# Read trimming & QC
# =============================================================================
def _fastp_input(wildcards):
    """Return fastp input dict based on library type."""
    if IS_PAIRED:
        return {"r1": f"{RAW_DIR}/{wildcards.sample}{R1_SUFFIX}",
                "r2": f"{RAW_DIR}/{wildcards.sample}{R2_SUFFIX}"}
    return {"r1": f"{RAW_DIR}/{wildcards.sample}{SE_SUFFIX}"}

def _fastp_output(wildcards):
    """Return fastp output dict based on library type."""
    d = {"html": f"{TRIM_DIR}/{wildcards.sample}.fastp.html",
         "json": f"{TRIM_DIR}/{wildcards.sample}.fastp.json"}
    if IS_PAIRED:
        d["r1"] = f"{TRIM_DIR}/{wildcards.sample}_1.trim.fq.gz"
        d["r2"] = f"{TRIM_DIR}/{wildcards.sample}_2.trim.fq.gz"
    else:
        d["r1"] = f"{TRIM_DIR}/{wildcards.sample}.trim.fq.gz"
    return d

rule fastp:
    input:
        unpack(_fastp_input)
    output:
        unpack(_fastp_output)
    conda: "envs/ge_analysis.yaml"
    threads: THREADS_FASTP
    run:
        if IS_PAIRED:
            shell(
                "fastp -i {input.r1} -I {input.r2} "
                "-o {output.r1} -O {output.r2} "
                "--html {output.html} --json {output.json} "
                "--thread {threads}"
            )
        else:
            shell(
                "fastp -i {input.r1} -o {output.r1} "
                "--html {output.html} --json {output.json} "
                "--thread {threads}"
            )

def _fastqc_input(wildcards):
    """Return fastqc input dict based on library type."""
    if IS_PAIRED:
        return {"r1": f"{TRIM_DIR}/{wildcards.sample}_1.trim.fq.gz",
                "r2": f"{TRIM_DIR}/{wildcards.sample}_2.trim.fq.gz"}
    return {"r1": f"{TRIM_DIR}/{wildcards.sample}.trim.fq.gz"}

def _fastqc_output(wildcards):
    """Return fastqc output dict based on library type."""
    if IS_PAIRED:
        return {"r1_html": f"{FASTQC_DIR}/{wildcards.sample}_1.trim_fastqc.html",
                "r2_html": f"{FASTQC_DIR}/{wildcards.sample}_2.trim_fastqc.html",
                "r1_zip":  f"{FASTQC_DIR}/{wildcards.sample}_1.trim_fastqc.zip",
                "r2_zip":  f"{FASTQC_DIR}/{wildcards.sample}_2.trim_fastqc.zip"}
    return {"html": f"{FASTQC_DIR}/{wildcards.sample}.trim_fastqc.html",
            "zip":  f"{FASTQC_DIR}/{wildcards.sample}.trim_fastqc.zip"}

rule fastqc:
    input:
        unpack(_fastqc_input)
    output:
        unpack(_fastqc_output)
    conda: "envs/ge_analysis.yaml"
    threads: THREADS_FASTQC
    run:
        if IS_PAIRED:
            shell(
                "fastqc {input.r1} {input.r2} "
                "--outdir {FASTQC_DIR} --threads {threads}"
            )
        else:
            shell(
                "fastqc {input.r1} "
                "--outdir {FASTQC_DIR} --threads {threads}"
            )

def _multiqc_input():
    """Return multiqc input list conditioned on library type."""
    if IS_PAIRED:
        return expand(f"{FASTQC_DIR}/{{sample}}_1.trim_fastqc.html", sample=SAMPLES)
    return expand(f"{FASTQC_DIR}/{{sample}}.trim_fastqc.html", sample=SAMPLES)

rule multiqc:
    input:
        unpack(_multiqc_input)
    output:
        report = f"{MULTIQC_DIR}/multiqc_report.html"
    conda: "envs/ge_analysis.yaml"
    shell:
        """
        multiqc {FASTQC_DIR} {TRIM_DIR} -o {MULTIQC_DIR} -f
        """

# =============================================================================
# Alignment
# =============================================================================
def _star_align_input(wildcards):
    """Return star_align input dict based on library type."""
    d = {"index": INDEX_DIR}
    if IS_PAIRED:
        d["r1"] = f"{TRIM_DIR}/{wildcards.sample}_1.trim.fq.gz"
        d["r2"] = f"{TRIM_DIR}/{wildcards.sample}_2.trim.fq.gz"
    else:
        d["r1"] = f"{TRIM_DIR}/{wildcards.sample}.trim.fq.gz"
    return d

rule star_align:
    input:
        unpack(_star_align_input)
    output:
        bam    = f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam",
        counts = f"{ALIGN_DIR}/{{sample}}_ReadsPerGene.out.tab",
        log    = f"{ALIGN_DIR}/{{sample}}_Log.final.out"
    threads: THREADS_STAR_ALIGN
    resources:
        mem_mb = MEM_MB_STAR_ALIGN
    conda: "envs/ge_analysis.yaml"
    params:
        prefix = f"{ALIGN_DIR}/{{sample}}_",
        sort_ram = MEM_MB_STAR_ALIGN * 1000000,   # STAR expects bytes
        tmpdir = lambda wildcards: f"/tmp/STAR_{wildcards.sample}"
    run:
        reads = "{input.r1} {input.r2}" if IS_PAIRED else "{input.r1}"
        shell(
            "rm -rf {params.tmpdir} && "
            "STAR --runThreadN {threads} "
            "--genomeDir {input.index} "
            "--readFilesIn " + reads + " "
            "--readFilesCommand gunzip -c "
            "--outFileNamePrefix {params.prefix} "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMunmapped Within "
            "--quantMode GeneCounts "
            "--outSAMattributes Standard "
            "--limitBAMsortRAM {params.sort_ram} "
            "--outTmpDir {params.tmpdir} && "
            "rm -rf {params.tmpdir}"
        )

rule samtools_index:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam"
    output:
        bai = f"{ALIGN_DIR}/{{sample}}_BAMIndex.bai"
    conda: "envs/ge_analysis.yaml"
    shell:
        """
        samtools index {input.bam} {output.bai}
        """

# =============================================================================
# Strandedness inference (RSeQC)
# =============================================================================
rule make_genes_bed:
    input:
        COMBINED_GTF_BED
    output:
        GENES_BED
    conda: "envs/ge_analysis.yaml"
    shell:
        """
        python scripts/gtf2bed.py {input} > {output}
        """

rule infer_strandedness:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam",
        bai = f"{ALIGN_DIR}/{{sample}}_BAMIndex.bai",
        bed = GENES_BED
    output:
        txt = f"{RSEQC_DIR}/{{sample}}_infer_experiment.txt"
    conda: "envs/ge_analysis.yaml"
    shell:
        r"""
        mkdir -p {RSEQC_DIR}
        infer_experiment.py -i {input.bam} -r {input.bed} > {output.txt}
        """

rule strandedness_summary:
    input:
        reports = expand(f"{RSEQC_DIR}/{{sample}}_infer_experiment.txt", sample=SAMPLES)
    output:
        tsv = f"{RSEQC_DIR}/strandedness_summary.tsv"
    conda: "envs/ge_analysis.yaml"
    run:
        import os

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        with open(output.tsv, 'w') as out_f:
            out_f.write("sample\tline\n")
            for f_path in input.reports:
                filename = os.path.basename(f_path)
                sample_name = filename.replace("_infer_experiment.txt", "")
                with open(f_path, 'r') as in_f:
                    for line in in_f:
                        line = line.strip()
                        if line.startswith("This is") or line.startswith("Fraction of reads"):
                            out_f.write(f"{sample_name}\t{line}\n")

# =============================================================================
# Count matrix generation
# =============================================================================
rule count_matrix:
    input:
        counts       = expand(f"{ALIGN_DIR}/{{sample}}_ReadsPerGene.out.tab", sample=SAMPLES),
        strandedness = f"{RSEQC_DIR}/strandedness_summary.tsv"
    output:
        matrix = f"{COUNTS_DIR}/gene_counts_matrix.tsv"
    conda: "envs/ge_analysis.yaml"
    run:
        import pandas as pd
        import os
        import re

        def parse_strandedness_summary(strandedness_file):
            """
            Parse strandedness summary and return a dict mapping sample -> column_name
            ('unstranded', 'first_strand', or 'second_strand').

            Logic for both PE and SE data:
              - If fractions are available and one strand clearly dominates
                (difference > 0.2), choose that strand.
              - Otherwise default to 'unstranded'.
            """
            sample_info = {}

            def decide_strand(first_frac, second_frac):
                """Decide column based on strand fractions."""
                if first_frac is not None and second_frac is not None:
                    diff = abs(first_frac - second_frac)
                    if diff > 0.2:
                        return "first_strand" if first_frac > second_frac else "second_strand"
                return "unstranded"

            with open(strandedness_file, 'r') as f:
                current_sample = None
                first_strand_frac = None
                second_strand_frac = None

                for line in f:
                    if line.startswith("sample"):
                        continue

                    parts = line.strip().split("\t", 1)
                    if len(parts) < 2:
                        continue

                    sample, content = parts

                    if current_sample is not None and sample != current_sample:
                        sample_info[current_sample] = decide_strand(first_strand_frac, second_strand_frac)
                        first_strand_frac = None
                        second_strand_frac = None

                    current_sample = sample

                    if 'Fraction of reads explained by "1++,1--,2+-,2-+"' in content:
                        match = re.search(r':\s*([\d.]+)', content)
                        if match:
                            first_strand_frac = float(match.group(1))

                    if 'Fraction of reads explained by "1+-,1-+,2++,2--"' in content:
                        match = re.search(r':\s*([\d.]+)', content)
                        if match:
                            second_strand_frac = float(match.group(1))

                if current_sample is not None:
                    sample_info[current_sample] = decide_strand(first_strand_frac, second_strand_frac)

            return sample_info

        column_mapping = parse_strandedness_summary(input.strandedness)
        counts_dict = {}

        for f in input.counts:
            sample = os.path.basename(f).replace("_ReadsPerGene.out.tab", "")
            df = pd.read_csv(f, sep="\t", header=None, index_col=0,
                             names=["gene_id", "unstranded", "first_strand", "second_strand"],
                             skiprows=4)
            column_name = column_mapping.get(sample, "unstranded")
            counts_dict[sample] = df[column_name]
            print(f"Sample {sample}: using {column_name} column")

        matrix_df = pd.DataFrame(counts_dict)
        matrix_df.index.name = "gene_id"
        os.makedirs(os.path.dirname(output.matrix), exist_ok=True)
        matrix_df.to_csv(output.matrix, sep="\t")

# =============================================================================
# Separate human and viral counts
# =============================================================================
rule separate_viral_counts:
    input:
        combined_matrix = f"{COUNTS_DIR}/gene_counts_matrix.tsv",
        pathogen_gtf    = PATHOGEN_GTF
    output:
        human_matrix = f"{COUNTS_DIR}/gene_counts_matrix_human.tsv",
        viral_matrix = f"{COUNTS_DIR}/gene_counts_matrix_viral.tsv"
    conda: "envs/ge_analysis.yaml"
    run:
        import pandas as pd
        import gzip
        import os

        matrix_df = pd.read_csv(input.combined_matrix, sep="\t", index_col=0)

        # Extract pathogen gene IDs from the GTF file (handle gzipped or plain text)
        viral_gene_ids = set()
        _is_gzipped = input.pathogen_gtf.endswith('.gz')
        if not _is_gzipped:
            try:
                with open(input.pathogen_gtf, 'rb') as _check:
                    _is_gzipped = _check.read(2) == b'\x1f\x8b'
            except (OSError, IOError):
                pass
        _open = gzip.open if _is_gzipped else open
        with _open(input.pathogen_gtf, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    attributes = parts[8]
                    if 'gene_id' in attributes:
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if attr.startswith('gene_id'):
                                gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                                viral_gene_ids.add(gene_id)

        viral_indices = []
        human_indices = []

        for gene_id in matrix_df.index:
            is_viral = False
            if gene_id in viral_gene_ids:
                is_viral = True
            else:
                gid_str = str(gene_id)
                for prefix in config.get("pathogen_gene_prefixes", []):
                    if prefix and gid_str.startswith(prefix):
                        is_viral = True
                        break

            if is_viral:
                viral_indices.append(gene_id)
            else:
                human_indices.append(gene_id)

        human_df = matrix_df.loc[human_indices] if human_indices else pd.DataFrame()
        viral_df = matrix_df.loc[viral_indices] if viral_indices else pd.DataFrame()

        os.makedirs(os.path.dirname(output.human_matrix), exist_ok=True)
        human_df.to_csv(output.human_matrix, sep="\t")
        viral_df.to_csv(output.viral_matrix, sep="\t")

        print(f"Separated counts: {len(human_indices)} human genes, {len(viral_indices)} viral genes")

# =============================================================================
# Differential expression – DESeq2
# =============================================================================
rule deseq2_analysis:
    input:
        counts  = f"{COUNTS_DIR}/gene_counts_matrix.tsv",
        samples = SAMPLES_TSV
    output:
        results          = f"{RESULTS_DIR}/deg_results.tsv",
        results_filtered = f"{RESULTS_DIR}/deg_results_significant.tsv",
        normalized_counts = f"{RESULTS_DIR}/normalized_counts.tsv",
        pca      = f"{RESULTS_DIR}/pca_plot.png",
        distance = f"{RESULTS_DIR}/distance_matrix.png",
        volcano  = f"{RESULTS_DIR}/volcano_plot.png",
        heatmap  = f"{RESULTS_DIR}/heatmap.png"
    log:
        f"{LOGS_DIR}/deseq2.log"
    conda: "envs/ge_analysis.yaml"
    script:
        "scripts/deseq2_analysis.R"

rule deseq2_analysis_human:
    input:
        counts  = f"{COUNTS_DIR}/gene_counts_matrix_human.tsv",
        samples = SAMPLES_TSV
    output:
        results          = f"{RESULTS_DIR}/deg_results_human.tsv",
        results_filtered = f"{RESULTS_DIR}/deg_results_significant_human.tsv",
        normalized_counts = f"{RESULTS_DIR}/normalized_counts_human.tsv",
        pca      = f"{RESULTS_DIR}/pca_plot_human.png",
        distance = f"{RESULTS_DIR}/distance_matrix_human.png",
        volcano  = f"{RESULTS_DIR}/volcano_plot_human.png",
        heatmap  = f"{RESULTS_DIR}/heatmap_human.png"
    log:
        f"{LOGS_DIR}/deseq2_human.log"
    conda: "envs/ge_analysis.yaml"
    script:
        "scripts/deseq2_analysis.R"

rule deseq2_analysis_viral:
    input:
        counts  = f"{COUNTS_DIR}/gene_counts_matrix_viral.tsv",
        samples = SAMPLES_TSV
    output:
        results          = f"{RESULTS_DIR}/deg_results_viral.tsv",
        results_filtered = f"{RESULTS_DIR}/deg_results_significant_viral.tsv",
        normalized_counts = f"{RESULTS_DIR}/normalized_counts_viral.tsv",
        pca      = f"{RESULTS_DIR}/pca_plot_viral.png",
        distance = f"{RESULTS_DIR}/distance_matrix_viral.png",
        volcano  = f"{RESULTS_DIR}/volcano_plot_viral.png",
        heatmap  = f"{RESULTS_DIR}/heatmap_viral.png"
    log:
        f"{LOGS_DIR}/deseq2_viral.log"
    conda: "envs/ge_analysis.yaml"
    script:
        "scripts/deseq2_analysis_viral.R"

# =============================================================================
# Pathway / GO enrichment (clusterProfiler)
# =============================================================================
rule enrichment_analysis:
    input:
        f"{RESULTS_DIR}/deg_results_human.tsv"
    output:
        f"{RESULTS_DIR}/enrichment/enrich_KEGG_GSEA.tsv",
        f"{RESULTS_DIR}/enrichment/enrich_GOBP_GSEA.tsv",
        f"{RESULTS_DIR}/enrichment/enrich_KEGG_GSEA_dotplot.png",
        f"{RESULTS_DIR}/enrichment/enrich_GOBP_GSEA_dotplot.png",
        f"{RESULTS_DIR}/enrichment/enrich_KEGG_ORA.tsv",
        f"{RESULTS_DIR}/enrichment/enrich_GOBP_ORA.tsv",
        f"{RESULTS_DIR}/enrichment/enrich_KEGG_ORA_dotplot.png",
        f"{RESULTS_DIR}/enrichment/enrich_GOBP_ORA_dotplot.png"
    log:
        f"{LOGS_DIR}/enrichment.log"
    conda: "envs/enrichment.yaml"
    params:
        pvalue   = config.get("enrichment_pvalue", 0.2),
        qvalue   = config.get("enrichment_qvalue", 0.2),
        organism = config.get("kegg_organism", "hsa"),
        org_db   = config.get("org_db", "org.Hs.eg.db"),
        id_type  = config.get("gene_id_type", "ENSEMBL")
    shell:
        "Rscript scripts/enrichment_analysis_standalone.R {params.pvalue} {params.qvalue} {params.organism} {params.org_db} {params.id_type} > {log} 2>&1"

# =============================================================================
# Summary HTML report
# =============================================================================
rule generate_report:
    input:
        # Human DESeq2
        deg_human      = f"{RESULTS_DIR}/deg_results_human.tsv",
        deg_sig_human  = f"{RESULTS_DIR}/deg_results_significant_human.tsv",
        norm_human     = f"{RESULTS_DIR}/normalized_counts_human.tsv",
        pca_human      = f"{RESULTS_DIR}/pca_plot_human.png",
        volcano_human  = f"{RESULTS_DIR}/volcano_plot_human.png",
        heatmap_human  = f"{RESULTS_DIR}/heatmap_human.png",
        distance_human = f"{RESULTS_DIR}/distance_matrix_human.png",
        # Viral DESeq2
        deg_viral      = f"{RESULTS_DIR}/deg_results_viral.tsv",
        deg_sig_viral  = f"{RESULTS_DIR}/deg_results_significant_viral.tsv",
        norm_viral     = f"{RESULTS_DIR}/normalized_counts_viral.tsv",
        pca_viral      = f"{RESULTS_DIR}/pca_plot_viral.png",
        volcano_viral  = f"{RESULTS_DIR}/volcano_plot_viral.png",
        heatmap_viral  = f"{RESULTS_DIR}/heatmap_viral.png",
        distance_viral = f"{RESULTS_DIR}/distance_matrix_viral.png",
        # Combined
        deg_combined   = f"{RESULTS_DIR}/deg_results.tsv",
        norm_combined  = f"{RESULTS_DIR}/normalized_counts.tsv",
        # Enrichment
        enrich_kegg        = f"{RESULTS_DIR}/enrichment/enrich_KEGG_GSEA.tsv",
        enrich_gobp        = f"{RESULTS_DIR}/enrichment/enrich_GOBP_GSEA.tsv",
        enrich_kegg_dot    = f"{RESULTS_DIR}/enrichment/enrich_KEGG_GSEA_dotplot.png",
        enrich_gobp_dot    = f"{RESULTS_DIR}/enrichment/enrich_GOBP_GSEA_dotplot.png",
        enrich_kegg_ora    = f"{RESULTS_DIR}/enrichment/enrich_KEGG_ORA.tsv",
        enrich_gobp_ora    = f"{RESULTS_DIR}/enrichment/enrich_GOBP_ORA.tsv",
        enrich_kegg_ora_dot= f"{RESULTS_DIR}/enrichment/enrich_KEGG_ORA_dotplot.png",
        enrich_gobp_ora_dot= f"{RESULTS_DIR}/enrichment/enrich_GOBP_ORA_dotplot.png",
        # Counts
        count_matrix  = f"{COUNTS_DIR}/gene_counts_matrix.tsv",
        count_human   = f"{COUNTS_DIR}/gene_counts_matrix_human.tsv",
        count_viral   = f"{COUNTS_DIR}/gene_counts_matrix_viral.tsv",
        # QC
        multiqc       = f"{MULTIQC_DIR}/multiqc_report.html",
        strandedness  = f"{RSEQC_DIR}/strandedness_summary.tsv",
        # Versions
        versions      = f"{RESULTS_DIR}/versions.txt"
    output:
        report = f"{RESULTS_DIR}/report.html"
    log:
        f"{LOGS_DIR}/report.log"
    conda: "envs/ge_analysis.yaml"
    script:
        "scripts/generate_report.R"

# =============================================================================
# Version log for publication reproducibility
# =============================================================================
rule version_log:
    """Write tool and package versions to results/versions.txt."""
    output:
        f"{RESULTS_DIR}/versions.txt"
    log:
        f"{LOGS_DIR}/versions.log"
    conda: "envs/ge_analysis.yaml"
    shell:
        r"""
        echo "VirusHostRNA — Tool Versions" > {output}
        echo "=============================" >> {output}
        echo "" >> {output}
        echo "--- Core Tools ---" >> {output}
        echo "STAR:"    $$(STAR --version 2>&1 | head -1) >> {output}
        echo "fastp:"   $$(fastp --version 2>&1 | head -1) >> {output}
        echo "samtools:" $$(samtools --version 2>&1 | head -1) >> {output}
        echo "fastqc:"  $$(fastqc --version 2>&1)            >> {output}
        echo "multiqc:" $$(multiqc --version 2>&1)          >> {output}
        echo "Python:"  $$(python --version 2>&1)           >> {output}
        echo "R:"       $$(R --version 2>&1 | head -1)      >> {output}
        echo "" >> {output}
        echo "--- R Packages ---" >> {output}
        Rscript -e 'cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")' >> {output}
        Rscript -e 'cat("clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")' >> {output}
        Rscript -e 'cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")' >> {output}
        Rscript -e 'cat("biomaRt:", as.character(packageVersion("biomaRt")), "\n")' >> {output}
        Rscript -e 'cat("pheatmap:", as.character(packageVersion("pheatmap")), "\n")' >> {output}
        Rscript -e 'cat("ggrepel:", as.character(packageVersion("ggrepel")), "\n")' >> {output}
        echo "" >> {output}
        echo "--- Pipeline ---" >> {output}
        echo "Snakemake:" $$(snakemake --version 2>&1) >> {output}
        echo "" >> {output}
        echo "Generated:" $$(date) >> {output}
        """ 2> {log}

# =============================================================================
# Convenience: mark alignments as complete without re-running STAR
# Run manually: snakemake touch_complete --cores 1
# =============================================================================
rule touch_complete:
    """Skip re-alignment — mark existing BAMs and counts as up-to-date."""
    output:
        touch(expand(f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)),
        touch(expand(f"{ALIGN_DIR}/{{sample}}_ReadsPerGene.out.tab", sample=SAMPLES))
    shell:
        "echo 'BAMs and gene counts marked as complete.'"

# =============================================================================
# HiPathia pathway analysis (commented out – enable if needed)
# =============================================================================
#rule hipathia_analysis:
#    input:
#        samples = SAMPLES_TSV
#    output:
#        pathway_activity = f"{RESULTS_DIR}/hipathia_pathway_activity.tsv",
#        pathway_heatmap  = f"{RESULTS_DIR}/hipathia_pathway_heatmap.png",
#        report_html      = f"{RESULTS_DIR}/hipathia_report/index.html"
#    log:
#        f"{LOGS_DIR}/hipathia.log"
#    conda: "envs/hipathia.yaml"
#    script:
#        "scripts/hipathia_analysis.R"
