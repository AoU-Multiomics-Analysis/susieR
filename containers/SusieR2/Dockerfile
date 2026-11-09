FROM mambaorg/micromamba:1.5.3@sha256:bb9f5c6a41682f0f77a80784d82a51f1947224d3c0e256e74c32c77c53d3f63e

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib
ENV SUSIER_SCRIPTS_PATH=/opt/r/scripts
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 VROOM_THREADS=1
ENV SUSIER_SUSIER_SOURCE_SHA=ef213feed2cb82419677661a8c986e1504df2c73
ENV SUSIER_EQTLUTILS_SOURCE_SHA=26242562a4e244334fd9691d03bc1ef4d2d6c1d9

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/SusieR2/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --strict-channel-priority --file /tmp/environment.yml && \
    micromamba list --name base --json > /opt/conda/susier-v2-conda-packages.json && \
    micromamba clean --all --yes

# Check native headers before downloading or compiling the R source packages.
RUN test -r /opt/conda/include/zlib.h && \
    test -r /opt/conda/include/zstd.h && \
    test -r /opt/conda/include/libdeflate.h && \
    pkg-config --modversion zlib libzstd libdeflate

# Install only these checked archives. Conda supplies all their dependencies.
# Upstream tag 0.16.6 points to the fixed susieR commit below.
RUN Rscript -e ' \
    options(timeout = 600); \
    sources <- list( \
      list("eQTLUtils_0.0.0.9000.tar.gz", "https://codeload.github.com/kauralasoo/eQTLUtils/tar.gz/26242562a4e244334fd9691d03bc1ef4d2d6c1d9", "bd4e766adb511960c75b89fe4ca926ff1215ff75e1ee720dd5172793eb3f464a"), \
      list("susieR_0.16.6.tar.gz", "https://codeload.github.com/stephenslab/susieR/tar.gz/ef213feed2cb82419677661a8c986e1504df2c73", "f36dd4045a448c0b07a1ea8429dc531d573e31fbe822bfe61842244ddd5ddd88"), \
      list("pgenlibr_0.7.0.tar.gz", c("https://cloud.r-project.org/src/contrib/pgenlibr_0.7.0.tar.gz", "https://cloud.r-project.org/src/contrib/Archive/pgenlibr/pgenlibr_0.7.0.tar.gz"), "b9ee17643b41d21053a65e2349e4f1e737ca07caa404628bed720510554a3858") \
    ); \
    for (source in sources) { \
      archive <- file.path(tempdir(), source[[1]]); \
      message("Installing checked source: ", source[[1]]); \
      downloaded <- FALSE; \
      for (url in source[[2]]) { \
        status <- tryCatch(download.file(url, archive, mode = "wb", method = "libcurl"), error = function(error) 1L); \
        if (isTRUE(status == 0L)) { downloaded <- TRUE; break }; \
      }; \
      stopifnot(downloaded); \
      stopifnot(identical(digest::digest(file = archive, algo = "sha256"), source[[3]])); \
      install.packages(archive, repos = NULL, type = "source", dependencies = FALSE); \
      unlink(archive); \
    }; \
    stopifnot(as.character(packageVersion("susieR")) == "0.16.6", \
              as.character(packageVersion("pgenlibr")) == "0.7.0", \
              as.character(packageVersion("eQTLUtils")) == "0.0.0.9000"); \
    invisible(lapply(c("tidyverse", "data.table", "optparse", "arrow", "devtools", \
                       "Rfast", "matrixStats", "digest", "jsonlite", "caret", \
                       "SummarizedExperiment", "Rsamtools", "eQTLUtils", "susieR", \
                       "pgenlibr", "R.utils"), function(package) stopifnot(requireNamespace(package, quietly = TRUE))))'

USER root
RUN mkdir -p /opt/r/lib /opt/r/scripts
COPY R/utils/ImportFunctions.R \
     R/utils/InitFunctions.R \
     R/utils/MemoryHelpers.R \
     R/utils/OptParser.R \
     R/utils/SusiePerformance.R \
     R/utils/SusieFunctions.R \
     R/utils/PgenFunctions.R \
     R/utils/SusieOutputFunctions.R \
     R/utils/PgenGeneFunctions.R \
     R/utils/PgenBatchFunctions.R \
     R/utils/PgenQcFunctions.R /opt/r/lib/
COPY R/scripts/plan_pgen_susie.R \
     R/scripts/filter_pgen_susie.R \
     R/scripts/run_pgen_susie_batch.R \
     R/scripts/aggregate_pgen_susie.R \
     R/scripts/susie.R /opt/r/scripts/
USER $MAMBA_USER

RUN plink2 --version && bgzip --version && tabix --version
CMD ["bash"]
