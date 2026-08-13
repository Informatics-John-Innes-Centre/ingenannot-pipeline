# Ingenannot Pipeline

## Overview 
A nextflow pipeline for genome annotation. The pipeline requires:
* A reference genome.
* Paired-end Illumina RNA-Seq data.
* PacBio Iso-Seq data.
* A protein database.

It creates annotations using 4 core annotators: [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Helixer](https://github.com/usadellab/Helixer), [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) and [Tiberius](https://github.com/Gaius-Augustus/Tiberius). It then uses [InGenAnnot](https://forge.inrae.fr/bioger/ingenannot) to compare and combine them to create an optimal single result annotation. The pipeline is designed such that it can process multiple seperate genomes with their own corresponding data at once.

## Preparing Input Data

A recommended way to arrange your input data is to create 3 folders:

```bash
input
├── genomes
├── isoseq
├── rnaseq
└── protein_database.fa
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
> The file endings of your input data have to be exact: 
> * Reference genomes **must** be `.fasta` files.    
> * RNA-Seq data **must** end in either `-r1.fastq.gz` or `-r2.fastq.gz`.    
> * The Iso-Seq data **must** be `.flnc.cram`.    

Create a file called `input.csv` in the following format:
```
genome_prefix,illumina_prefix,isoseq_prefix
genomeA,genomeA,genomeA
genomeB,genomeB,genomeB
```

Each row is a seperate genome. It is worth noting that the prefixes are unique to each type of data and therefore, **do not have to be the same**.

## Running
Clone the repository and move your inputs (input data and CSV file) into it:

```bash
git clone github.com/Informatics-John-Innes-Centre/ingenannot-pipeline
cd ingenannot-pipeline
mv /path/to/my/inputs .
mv /path/to/my/input.csv .
```

To install required dependencies you will need to first install the [Pixi](https://github.com/prefix-dev/pixi/) package manager. This pipeline requires **nextflow** and **apptainer**. You can opt to use your own installation of apptainer like so (recommended for **HPC environments): 
```bash
pixi install
```
Or you can include apptainer in the pixi environment (recommended for **local** execution):
```bash
pixi install -e apptainer
```

Build all the containers required for the pipeline with the `containers.sh` script (this will take a while):
```bash
./containers.sh
```

### Local

```bash
pixi run nextflow run main.nf -profile local \
    --frozDir "/path/to/your/genomes" \
    --accessionFile "/path/to/your/input.csv" \
    --proteinDatabase "/path/to/your/protein_database.fa" \
    --fastqDirectory "/path/to/your/rnaseq" \
    --isoseqDirectory "/path/to/your/isoseq" \
    --tiberiusModel "tiberius_model" \
    --helixerLineage "helixer_lineage" \
    --annevoModel "annevo_model" \
    --annevoLineage "annevo_lineage" \
```
### Slurm

```bash
pixi run nextflow run main.nf -profile hpc \
    --frozDir "/path/to/your/genomes" \
    --accessionFile "/path/to/your/input.csv" \
    --proteinDatabase "/path/to/your/protein_database.fa" \
    --fastqDirectory "/path/to/your/rnaseq" \
    --isoseqDirectory "/path/to/your/isoseq" \
    --tiberiusModel "tiberius_model" \
    --helixerLineage "helixer_lineage" \
    --annevoModel "annevo_model" \
    --annevoLineage "annevo_lineage" \
    --slurm_queue your-queue \
    --slurm_gpu_queue your-queue \
```

`slurm_gpu_queue` is the queue that all GPU accelerated tasks will run on, all others will run on `slurm_queue`.
