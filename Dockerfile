FROM mambaorg/micromamba:1.5.10

USER root
WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

USER $MAMBA_USER
RUN micromamba install -y -n base -f /tmp/environment.yml \
        -c conda-forge -c bioconda diamond curl && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ /app/scripts/
COPY --chown=$MAMBA_USER:$MAMBA_USER webapp/ /app/webapp/
COPY --chown=$MAMBA_USER:$MAMBA_USER batch_dehost.sh /app/batch_dehost.sh

RUN chmod +x /app/scripts/bootstrap_offline_data.sh /app/batch_dehost.sh

ENV PATH="/opt/conda/bin:$PATH" \
    CDC_INPUT_DIR="/data/input" \
    CDC_OUTPUT_DIR="/data/output" \
    CDC_DB_DIR="/data/databases" \
    CDC_REF_DIR="/data/references" \
    CDC_PORT="8080" \
    CDC_AUTO_DOWNLOAD_K2="false"

EXPOSE 8080

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/app/scripts/bootstrap_offline_data.sh"]
CMD ["python3", "/app/webapp/server.py"]

