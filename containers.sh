#! /bin/bash

# Pull biocontainers
mkdir -p singularity-cache
singularity pull --dir singularity-cache docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3
mv bedtools_2.31.1--h13024bc_3.sif quay.io-biocontainers-bedtools-2.31.1--h13024bc_3.img

# Build custom containers
sudo singularity build singularity-cache/helixer.sif containers/helixer.def
sudo singularity build singularity-cache/annevo.sif containers/annevo.def
sudo singularity build singularity-cache/tiberius.sif containers/tiberius.def
sudo singularity build singularity-cache/miniprot.sif containers/miniprot.def
sudo singularity build singularity-cache/cutadapt.sif containers/cutadapt.def
sudo singularity build singularity-cache/star.sif containers/star.def
sudo singularity build singularity-cache/braker3.sif containers/braker3.def
sudo singularity build singularity-cache/stringtie.sif containers/stringtie.def
sudo singularity build singularity-cache/ingenannot.sif containers/ingenannot.def
sudo singularity build singularity-cache/agat.sif containers/agat.def
sudo singularity build singularity-cache/samtools.sif containers/samtools.def

