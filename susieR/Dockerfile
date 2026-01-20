## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
## Create a new environment and install packages
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    conda-forge::r-tidyverse \
    conda-forge::datatable \
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
    conda-forge::r-qs2


RUN R -e "remotes::install_github('kauralasoo/eQTLUtils')"

COPY susieR/*  . 
COPY susieR/utils/*  /opt/r/lib/ 


## Default command to run when the container starts
CMD ["bash"]
