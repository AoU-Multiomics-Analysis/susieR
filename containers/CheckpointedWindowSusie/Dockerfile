FROM ghcr.io/aou-multiomics-analysis/susier@sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764

ENV CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST=sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764
ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib

USER root
COPY R/utils/CheckpointedWindowSusieFunctions.R R/utils/CheckpointStore.R /opt/r/lib/
COPY R/scripts/run_checkpointed_window_susie.R /opt/r/scripts/
RUN ln -sf /opt/r/scripts/run_checkpointed_window_susie.R /tmp/run_checkpointed_window_susie.R
USER $MAMBA_USER

CMD ["bash"]
