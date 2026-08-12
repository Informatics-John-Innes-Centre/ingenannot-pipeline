#! /bin/bash

# Pull biocontainers
mkdir -p apptainer-cache

# Build custom containers
pixi run apptainer build apptainer-cache/bedtools.sif containers/bedtools.def
pixi run apptainer build apptainer-cache/helixer.sif containers/helixer.def
pixi run apptainer build apptainer-cache/annevo.sif containers/annevo.def
pixi run apptainer build apptainer-cache/tiberius-1.1.8.sif containers/tiberius-1.1.8.def
pixi run apptainer build apptainer-cache/tiberius-2.0.6.sif containers/tiberius-2.0.6.def
pixi run apptainer build apptainer-cache/miniprot.sif containers/miniprot.def
pixi run apptainer build apptainer-cache/cutadapt.sif containers/cutadapt.def
pixi run apptainer build apptainer-cache/star.sif containers/star.def
pixi run apptainer build apptainer-cache/braker3.sif containers/braker3.def
pixi run apptainer build apptainer-cache/stringtie.sif containers/stringtie.def
pixi run apptainer build apptainer-cache/ingenannot.sif containers/ingenannot.def
pixi run apptainer build apptainer-cache/samtools.sif containers/samtools.def

