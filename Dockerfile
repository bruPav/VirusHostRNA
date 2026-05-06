# =============================================================================
# VirusHostRNA — reproducible container with all conda environments pre-built.
#
# Build:
#   docker build -t virushostrna .
#
# Run (mount your data + references):
#   docker run --rm -v $PWD:/data -v $PWD/refs:/opt/pipeline/refs \
#     virushostrna snakemake --use-conda --cores 8
#
# Singularity / Apptainer:
#   singularity build pipeline.sif docker-daemon://virushostrna:latest
#   singularity exec --bind $PWD:/data pipeline.sif \
#     snakemake --use-conda --cores 8 --conda-prefix /opt/conda/envs
# =============================================================================

FROM continuumio/miniforge3:latest

# Use mamba for faster environment solving
RUN conda install -n base -c conda-forge mamba -y && \
    conda clean -afy

# Copy environment definitions first (layer caching — rebuilds only when yamls change)
COPY envs/ /opt/pipeline/envs/

# Pre-build all conda environments
RUN mamba env create -f /opt/pipeline/envs/ge_analysis.yaml -q && \
    mamba env create -f /opt/pipeline/envs/enrichment.yaml -q && \
    mamba env create -f /opt/pipeline/envs/hipathia.yaml -q && \
    find /opt/conda/envs -name '*.tar.bz2' -delete && \
    conda clean -afy

# Install Snakemake into base environment
RUN mamba install -n base -c conda-forge -c bioconda snakemake -y && \
    conda clean -afy

# Copy pipeline code
COPY . /opt/pipeline/

WORKDIR /opt/pipeline

ENTRYPOINT ["snakemake"]
CMD ["--use-conda", "--conda-prefix", "/opt/conda/envs"]
