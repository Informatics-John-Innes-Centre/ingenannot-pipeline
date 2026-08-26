#! /bin/bash -e

# Pull biocontainers
mkdir -p apptainer-cache

# Build custom containers
pixi run apptainer pull apptainer-cache/bedtools.sif docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3
pixi run apptainer pull apptainer-cache/annevo.sif docker://quay.io/biocontainers/annevo:2.3.2--hdfd78af_0
pixi run apptainer pull apptainer-cache/tiberius.sif docker://larsgabriel23/tiberius:2.0.6
pixi run apptainer pull apptainer-cache/miniprot.sif docker://quay.io/biocontainers/miniprot:0.18--h577a1d6_0
pixi run apptainer pull apptainer-cache/cutadapt.sif docker://quay.io/biocontainers/cutadapt:5.2--py310h6813faf_2
pixi run apptainer pull apptainer-cache/star.sif docker://quay.io/biocontainers/star:2.7.11b--h5ca1c30_8
pixi run apptainer pull apptainer-cache/braker3.sif docker://teambraker/braker3:v3.0.8
pixi run apptainer pull apptainer-cache/stringtie.sif docker://quay.io/biocontainers/stringtie:3.0.3--h29c0135_0
pixi run apptainer pull apptainer-cache/samtools.sif docker://quay.io/biocontainers/samtools:1.23--h96c455f_0
pixi run apptainer pull apptainer-cache/isoseq.sif docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0
pixi run apptainer pull apptainer-cache/pbmm2.sif docker://quay.io/biocontainers/pbmm2:26.1.99--h9ee0642_0
pixi run apptainer build apptainer-cache/helixer.sif containers/helixer.def
pixi run apptainer build apptainer-cache/ingenannot.sif containers/ingenannot.def
