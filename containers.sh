#! /bin/bash

# Pull biocontainers
mkdir -p singularity-cache

# Build custom containers
sudo singularity build singularity-cache/bedtools.sif containers/bedtools.def
sudo singularity build singularity-cache/helixer.sif containers/helixer.def
sudo singularity build singularity-cache/annevo.sif containers/annevo.def
sudo singularity build singularity-cache/tiberius-1.1.8.sif containers/tiberius-1.1.8.def
sudo singularity build singularity-cache/tiberius-2.0.6.sif containers/tiberius-2.0.6.def
sudo singularity build singularity-cache/miniprot.sif containers/miniprot.def
sudo singularity build singularity-cache/cutadapt.sif containers/cutadapt.def
sudo singularity build singularity-cache/star.sif containers/star.def
sudo singularity build singularity-cache/braker3.sif containers/braker3.def
sudo singularity build singularity-cache/stringtie.sif containers/stringtie.def
sudo singularity build singularity-cache/ingenannot.sif containers/ingenannot.def
sudo singularity build singularity-cache/samtools.sif containers/samtools.def

