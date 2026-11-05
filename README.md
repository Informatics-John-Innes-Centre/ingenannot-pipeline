# Ingenannot Pipeline

## Overview 

A Nextflow pipeline for genome annotation. The pipeline requires:
* A genome.
* A repeat mask BED file for that genome, with only three columns. See [here](example/genomes/my_arabidopsis_all_repeats.bed) for an example.
* Paired-end Illumina RNA-Seq data.
* PacBio Iso-Seq data **(optional)**.
* A protein database.

The pipeline creates annotations using 4 core annotators: [ANNEVO](https://github.com/xjtu-omics/ANNEVO), [Helixer](https://github.com/usadellab/Helixer), [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) and [Tiberius](https://github.com/Gaius-Augustus/Tiberius). It then uses [InGenAnnot](https://forge.inrae.fr/bioger/ingenannot) to compare and combine them to create an optimal single result annotation. The pipeline is designed such that it can process multiple separate genomes with their own corresponding data at once, and can run easily on an **offline** HPC cluster.

## Preparing Input Data

The pipeline requires an `input.csv`, which looks like this:

```
genome_prefix,illumina_prefix,isoseq_prefix
genomeA,illuminaA,isoseqA
genomeB,illuminaB,isoseqB
```

> [!NOTE]
> If you do not want to provide PacBio Iso-Seq data for a specific genome, leave `isoseq_prefix` blank for that row.

Each row in this file is a separate genome to annotate. The pipeline determines which input data belongs to which genome by these prefixes. A recommended way to provide your input data is using this folder structure:

```bash
input
├── genomes
├── isoseq
├── rnaseq
├── protein_database.fa
└── input.csv
```

All data of each type goes in their corresponding folders. For example:

```bash
input
├── genomes
│   ├── genomeA.fasta
│   ├── genomeA_all_repeats.bed
│   ├── genomeB.fasta
│   └── genomeB_all_repeats.bed
├── isoseq
│   ├── isoseqA.flnc.cram
│   └── isoseqB.flnc.cram
└── rnaseq
    ├── illuminaA-sample1-r1.fastq.gz
    ├── illuminaA-sample1-r2.fastq.gz
    ├── illuminaA-sample2-r1.fastq.gz
    ├── illuminaA-sample2-r2.fastq.gz
    ├── illuminaB-sample1-r1.fastq.gz
    └── illuminaB-sample1-r2.fastq.gz
```

This example has two genomes: `genomeA` and `genomeB`. Each genome has a corresponding `.flnc.cram` file. `genomeA` has two rnaseq samples and `genomeB` has just one.

> [!IMPORTANT]
> The file endings of your input data have to be exact: 
> * Reference genomes **must** be a fasta file with name formatted like so: `{genome_prefix}.fasta`.
> * The repeat mask BED file **must** be of the format `{genome_prefix}_all_repeats.bed` and contain only 3 columns.
> * RNA-Seq data **must** be gzipped fastq files ending in either `-r1.fastq.gz` or `-r2.fastq.gz`.    
> * The Iso-Seq data **must** must be full length non chimeric reads stored as an unaligned cram file, with the name being formatted like so: `{isoseq_prefix}.flnc.cram`.    

## Running

This pipeline can be run uncontainerized, or via a container system such as Apptainer. Note that if you choose to run the pipeline via a container system, the jobs will only be able to run on the node you start it on. This reduces potential task parallelism on a multi-node HPC.

### Apptainer

Pull the container image like so:

```bash
apptainer pull docker://ghcr.io/informatics-john-innes-centre/ingenannot-pipeline:latest
```

Then you can run the pipeline:

```bash
apptainer run ingenannot-pipeline_latest.sif \
    --genomes "/path/to/your/genomes" \
    --csv "/path/to/your/input.csv" \
    --proteinDatabase "/path/to/your/protein_database.fa" \
    --rnaseq "/path/to/your/rnaseq" \
    --isoseq "/path/to/your/isoseq" \ # <-- optional
    --tiberiusModel "tiberius_model" \
    --helixerLineage "helixer_lineage" \
    --annevoModel "annevo_model" \
    --annevoLineage "annevo_lineage" \
```

### Uncontainerized

To run the pipeline uncontainerized, you will need **Apptainer** and **Nextflow**. First you need to populate the container cache for the pipeline:

```bash
./containers.sh
```

Then you can either run the pipeline **locally**:

```bash
nextflow run main.nf -profile local \
    --genomes "/path/to/your/genomes" \
    --csv "/path/to/your/input.csv" \
    --proteinDatabase "/path/to/your/protein_database.fa" \
    --rnaseq "/path/to/your/rnaseq" \
# optional
    --isoseq "/path/to/your/isoseq" \
    --tiberiusModel "tiberius_model" \
    --helixerLineage "helixer_lineage" \
    --annevoModel "annevo_model" \
    --annevoLineage "annevo_lineage" \
```

Or via a **SLURM** cluster:

```bash
nextflow run main.nf -profile hpc \
    --genomes "/path/to/your/genomes" \
    --csv "/path/to/your/input.csv" \
    --proteinDatabase "/path/to/your/protein_database.fa" \
    --rnaseq "/path/to/your/rnaseq" \
# optional
    --isoseq "/path/to/your/isoseq" \
    --tiberiusModel "tiberius_model" \
    --helixerLineage "helixer_lineage" \
    --annevoModel "annevo_model" \
    --annevoLineage "annevo_lineage" \
    --slurm_queue your-queue \
    --slurm_gpu_queue your-queue \
```

`slurm_gpu_queue` is the queue that all GPU accelerated tasks will run on, all others will run on `slurm_queue`.