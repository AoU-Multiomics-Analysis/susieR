## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3
#FROM ghcr.io/prefix-dev/pixi:latest

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin


## Create a new environment and install packages
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    conda-forge::r-tidyverse \
    conda-forge::r-data.table \
    conda-forge::r-optparse \
    conda-forge::r-arrow \
    bioconda::htslib \
    conda-forge::r-rnomni \
    conda-forge::gsutil \
    conda-forge::r-r.utils \
    conda-forge::r-janitor \
    conda-forge::r-magrittr \
    conda-forge::r-devtools \
    conda-forge::r-remotes \
    conda-forge::r-arrow \
    bioconda::bioconductor-genomicranges=1.58.0 \
    bioconda::bioconductor-rsamtools \
    bioconda::bioconductor-edger \
    bioconda::bioconductor-summarizedexperiment \
    conda-forge::r-susier \
    conda-forge::r-igraph \
    conda-forge::r-zoo \
    conda-forge::r-expm \
    conda-forge::r-matrixstats \
    conda-forge::r-rfast \
    conda-forge::r-assertthat \
    bioconda::bioconductor-pcatools \
    conda-forge::r-caret \
    conda-forge::r-stringfish \
    conda-forge::r-qs2 \
    conda-forge::r-lobstr

RUN R -e "remotes::install_github('kauralasoo/eQTLUtils')"
ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib
COPY R/utils/*  /opt/r/lib/
USER root
RUN mkdir -p /opt/r/scripts
COPY R/scripts/susie.R R/scripts/ComputeR2Susie.R R/scripts/PrepSusieCVPCs.R /opt/r/scripts/
RUN ln -s /opt/r/scripts/susie.R /susie.R && \
    ln -s /opt/r/scripts/ComputeR2Susie.R /ComputeR2Susie.R && \
    ln -s /opt/r/scripts/PrepSusieCVPCs.R /PrepSusieCVPCs.R && \
    ln -s /opt/r/scripts/susie.R /tmp/susie.R && \
    ln -s /opt/r/scripts/ComputeR2Susie.R /tmp/ComputeR2Susie.R && \
    ln -s /opt/r/scripts/PrepSusieCVPCs.R /tmp/PrepSusieCVPCs.R
USER $MAMBA_USER



## Default command to run when the container starts
CMD ["bash"]
