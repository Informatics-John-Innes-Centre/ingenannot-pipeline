# Ingenannot Pipeline

## Overview 
A nextflow pipeline for genome annotation. The pipeline requires:
* A reference genome.
* Paired-end Illumina RNA-Seq data.
* PacBio Iso-Seq data.

It creates annotations using 4 core annotators: [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Helixer](https://github.com/usadellab/Helixer), [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) and [Tiberius](https://github.com/Gaius-Augustus/Tiberius). It then uses [InGenAnnot](https://forge.inrae.fr/bioger/ingenannot) to compare and combine them to create an optimal single result annotation. The pipeline is designed such that it can process multiple seperate genomes with their own corresponding data at once.

## Preparing Input Data

A recommended way to arrange your input data is to create 3 folders:

```bash
input
├── genomes
├── isoseq
└── rnaseq
```

All the data of each type goes in their corresponding folders. The pipeline will determine which data belongs to which genome by a common prefix. For example:

```bash
input
├── genomes
│   ├── genomeA.fasta
│   └── genomeB.fasta
├── isoseq
│   ├── genomeA.flnc.cram
│   └── genomeB.flnc.cram
└── rnaseq
    ├── genomeA-sample1-r1.fastq.gz
    ├── genomeA-sample1-r2.fastq.gz
    ├── genomeA-sample2-r1.fastq.gz
    ├── genomeA-sample2-r2.fastq.gz
    ├── genomeB-sample1-r1.fastq.gz
    └── genomeB-sample1-r2.fastq.gz
```

This example has 2 genomes: `genomeA` and `genomeB`. Each genome has a corresponding `.flnc.cram` file. `genomeA` has two rnaseq samples and `genomeB` has just one.

> [!IMPORTANT]  
> Genomes **must** be `.fasta` files.    
> RNA-Seq data **must** end in either `-r1.fastq.gz` or `-r2.fastq.gz`.    
> The Iso-Seq data **must** be `.flnc.cram`.    

Create a file called `input.csv` in the following format:
```
genome_prefix,illumina_prefix,isoseq_prefix
genomeA,genomeA,genomeA
genomeB,genomeB,genomeB
```

Each row is a seperate genome. It is worth noting that the prefixes are unique to each type of data and therefore, **do not have to be the same**.
