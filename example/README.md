# Arabidopsis Example

This is an example dataset derived from the last 1,000,000 nucleotides of *Arabidopsis thalian's* chromosome Chr5, split into 8 synthetic contigs. The protein database is a subsampled version of a *Viridiplantae* OrthoDB plant protein database and there is a single pair of RNA-Seq paired-end fastq's. All of this data was derived from Braker 3's [example data](https://github.com/Gaius-Augustus/BRAKER#example-data), the only difference being the RNA-Seq bam had all non paired end reads filtered out and was then converted to a fastq pair.

## Runtime

On a HPC node with an AMD EPYC 7713P 64 core processor, 512GB of RAM, and **no GPU**, this example ran in 40 minutes. Peak memory usage was around 25 GB.