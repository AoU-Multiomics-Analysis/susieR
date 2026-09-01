FROM ghcr.io/aou-multiomics-analysis/susier@sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764

ENV CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST=sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764
ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib

RUN micromamba create -y -n checkpoint-gcloud -c conda-forge \
    google-cloud-sdk=581.0.0 && \
    micromamba clean --all --yes
ENV PATH=/opt/conda/envs/checkpoint-gcloud/bin:${PATH}

USER root
RUN ln -sf /opt/conda/envs/checkpoint-gcloud/bin/gsutil /opt/conda/bin/gsutil
COPY R/utils/CheckpointedWindowSusieFunctions.R R/utils/CheckpointStore.R /opt/r/lib/
COPY R/scripts/run_checkpointed_window_susie.R /opt/r/scripts/
RUN ln -sf /opt/r/scripts/run_checkpointed_window_susie.R /tmp/run_checkpointed_window_susie.R
USER $MAMBA_USER

CMD ["bash"]
