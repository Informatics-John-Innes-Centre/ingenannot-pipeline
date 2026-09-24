FROM debian:13.6

# Install base dependencies
RUN apt-get update && \
apt-get install -y wget \
ca-certificates \
locales-all \
procps && \
apt-get clean && \
rm -rf /var/lib/apt/lists/*

WORKDIR /opt/pipeline

# Copy stuff so we can install the environment
COPY pixi.lock pixi.toml .

# Install pixi environment
RUN wget -qO- https://pixi.sh/install.sh | PIXI_HOME=/usr/local sh && \
pixi config set mirrors '{"https://conda.anaconda.org/": ["https://prefix.dev/"]}' && \
pixi install -e apptainer

# Wrapper script
RUN printf '#!/bin/bash\n\nexec pixi run -e apptainer nextflow run /opt/pipeline/main.nf -c /opt/pipeline/nextflow.config -profile local "$@"' >> entrypoint.sh && \
chmod +x entrypoint.sh

# Copy nextflow script files
COPY main.nf nextflow.config .
ADD workflows workflows/

# Copy container images
# Intentionally copy one at a time to catch mistakes where one is missing from the folder
COPY apptainer-cache/samtools.sif \
apptainer-cache/helixer.sif \
apptainer-cache/ingenannot.sif \
apptainer-cache/tiberius.sif \
apptainer-cache/bedtools.sif \
apptainer-cache/annevo.sif \
apptainer-cache/miniprot.sif \
apptainer-cache/cutadapt.sif \
apptainer-cache/star.sif \
apptainer-cache/braker3.sif \
apptainer-cache/stringtie.sif \
apptainer-cache/samtools.sif \
apptainer-cache/isoseq.sif \
apptainer-cache/pbmm2.sif \
apptainer-cache/

# Environment stuff
ENV PIXI_PROJECT_MANIFEST=/opt/pipeline/pixi.toml
ENV NXF_OFFLINE='true'

# Run command
CMD [ "entrypoint.sh" ]