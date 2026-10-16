FROM debian:13.7
COPY main.nf /opt/pipeline/
ADD workflows /opt/pipeline/workflows
COPY nextflow.config /opt/pipeline/
ADD containers /opt/pipeline/containers
COPY containers.sh /opt/pipeline/
COPY pixi.lock /opt/pipeline/
COPY pixi.toml /opt/pipeline/
COPY apptainer-cache/samtools.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/helixer.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/ingenannot.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/tiberius.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/bedtools.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/annevo.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/miniprot.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/cutadapt.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/star.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/braker3.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/stringtie.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/samtools.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/isoseq.sif /opt/pipeline/apptainer-cache/
COPY apptainer-cache/pbmm2.sif /opt/pipeline/apptainer-cache/
ENV NXF_OFFLINE='true' \
PIXI_PROJECT_MANIFEST=/opt/pipeline/pixi.toml
WORKDIR /opt/pipeline
RUN apt-get update -y && \
apt-get install -y wget ca-certificates less which locales-all procps && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* && \
wget -qO- https://pixi.sh/install.sh | PIXI_HOME=/usr/local sh && \
pixi config set mirrors '{"https://conda.anaconda.org/": ["https://prefix.dev/"]}' && \
pixi install -e apptainer
ENTRYPOINT ["pixi", "run", "-e", "apptainer", "nextflow", "run", "/opt/pipeline/main.nf", "-c", "/opt/pipeline/nextflow.config", "-profile", "local"]