

process helixer {
    input:
    val genome_prefix
    path masked_fasta

    output:
    tuple val("helixer"), path("${genome_prefix}_helixer.gff")
    
    script:
    """
    Helixer.py\
        --downloaded-model-path /opt/helixer/models \
        --lineage land_plant \
        --fasta-path ${masked_fasta} \
        --species Pisum \
        --gff-output-path "${genome_prefix}_helixer.gff" \
        --peak-threshold 0.8 \
        --subsequence-length 64152 \
        --overlap-offset 32076 \
        --overlap-core-length 48114
    """
}

process annevo {
    input:
    val genome_prefix
    path masked_fasta

    output:
    tuple val("annevo"), path("${genome_prefix}_annevo.gff")
    
    script:
    """
    python /opt/ANNEVO/annotation.py \
        --genome ${masked_fasta} \
        --model_path /opt/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt \
        --output ${genome_prefix}_annevo.gff \
        --lineage Magnoliopsida \
        --threads ${task.cpus} \
        --num_workers ${task.cpus}
    """
}

process tiberius {
    input:
    val genome_prefix
    path masked_fasta

    output:
    tuple val("tiberius"), path("${genome_prefix}_tiberius.gtf")

    script:
    """
    python /opt/Tiberius/tiberius.py \
        --genome ${masked_fasta} \
        --model /opt/Tiberius/model_weights/eudicotyledons_weights \
        --out ${genome_prefix}_tiberius.gtf
    """
}

process miniprot  {
    input:
    val genome_prefix
    path masked_fasta

    output:
    path "${genome_prefix}_miniprot.gff3"

    script:
    """
    miniprot \
        -t ${task.cpus} \
        -d ${genome_prefix}_masked.mpi \
        ${masked_fasta}

    miniprot \
        -t ${task.cpus} \
        ${genome_prefix}_masked.mpi  \
        ${params.proteinDatabase} \
        --gff \
        > ${genome_prefix}_miniprot.gff3
    """
}

process star_index {
    input:
    val genome_prefix
    path masked_fasta

    output:
    path "${genome_prefix}_star_index/"

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --genomeDir ${genome_prefix}_star_index \
    --genomeFastaFiles ${masked_fasta} \
    --runThreadN ${task.cpus} \
    --limitGenomeGenerateRAM 90000000000 # Default is only 31G
    """
}

process trim_rnaseq_data {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(["trimmed-${reads[0]}", "trimmed-${reads[1]}"])

    script:
    """
    cutadapt \
        --cores ${task.cpus} \
        -m 75:75 \
        --nextseq-trim 20 \
        -a "ACTGTCTCTTATACACATCT" \
        -A "ACTGTCTCTTATACACATCT" \
        -o "trimmed-${reads[0]}" \
        -p "trimmed-${reads[1]}" \
        ${reads[0]} \
        ${reads[1]}
    """
}

process star_map_rnaseq_data_first_pass {
    input:
    tuple val(sample_id), path(reads)
    path star_index
    
    output:
    path("${sample_id}_SJ.out.tab")

    script:
    """
    STAR \
    --readFilesIn ${reads[0]} \
                    ${reads[1]} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${sample_id}_ \
    --runThreadN ${task.cpus} \
    --genomeDir ${star_index} \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 48000000000 \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNoverLmax 0.05 \
    --alignIntronMin 20 \
    --alignIntronMax 50000 \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --outFilterMultimapScoreRange 0 \
    --outFilterMatchNminOverLread 0.8 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --chimSegmentMin 50
    """
}

process star_map_rnaseq_data_second_pass {
    input:
    tuple val(sample_id), path(reads)
    path star_index
    path sj_files

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR \
    --readFilesIn ${reads[0]} \
                    ${reads[1]} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${sample_id}_ \
    --runThreadN ${task.cpus} \
    --genomeDir ${star_index} \
    --sjdbFileChrStartEnd  ${sj_files} \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 48000000000 \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNoverLmax 0.05 \
    --alignIntronMin 20 \
    --alignIntronMax 50000 \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --outFilterMultimapScoreRange 0 \
    --outFilterMatchNminOverLread 0.8 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --chimSegmentMin 50
    """
}


process braker3 {    
    input:
    val genome_prefix
    path masked_fasta
    path bam_files

    output:
    tuple val("braker3"), path("braker3.gff")

    script:
    """
    cp -r /opt/Augustus/config ./augustus_config
    export AUGUSTUS_CONFIG_PATH=\$PWD/augustus_config
    braker.pl \
    --prot_seq=${params.proteinDatabase} \
    --genome=${masked_fasta} \
    --bam=${bam_files.join(',')} \
    --threads=${task.cpus} \
    --species=${genome_prefix} \
    --workingdir=\$PWD \
    --gff3 \
    --verbosity=4
    """
}

process stringtie {
    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}.gtf"

    script:"""
    echo ${sample_id}
    echo ${bam}
    stringtie \
        -p ${task.cpus} \
        --nasc \
        ${bam} \
        > ${sample_id}.gtf \

    """
}

workflow annotate {
    take:
    masked_file
    genome_prefix
    illumina_prefix
    
    main:
    // Simple - annotators
    def helixer_annotation_ch = helixer(channel.value(genome_prefix), masked_file)
    def annevo_annotation_ch = annevo(channel.value(genome_prefix), masked_file)
    def tiberius_annotation_ch = tiberius(channel.value(genome_prefix), masked_file)
    
    // Miniprot alignment
    def miniprot_alignment_ch = miniprot(channel.value(genome_prefix), masked_file)
    // Star index used for star passes later
    def star_index = star_index(channel.value(genome_prefix), masked_file)
    // get all rna seq file pairs
    channel
        .fromFilePairs("${params.fastqDirectory}/${illumina_prefix}-*-r{1,2}.fastq.gz")
        .set { rnaseq_pairs }
    // trimming
    def trimmed_pairs = trim_rnaseq_data(rnaseq_pairs)
    // sj files from first pass
    def first_pass_mapped = star_map_rnaseq_data_first_pass(trimmed_pairs, star_index)
    def all_sj = first_pass_mapped.collect()
    // get bams from second pass
    def second_pass_bams = star_map_rnaseq_data_second_pass(trimmed_pairs, star_index, all_sj)
    def all_bams_ch = second_pass_bams
    .map { _sample_id, bam -> bam }
    def stringtie_transcripts_ch = stringtie(second_pass_bams)

    def braker_annotation_ch = braker3(channel.value(genome_prefix), masked_file, all_bams_ch.collect())    

    emit:
    helixer_annotation = helixer_annotation_ch
    annevo_annotation  = annevo_annotation_ch
    tiberius_annotation = tiberius_annotation_ch
    braker_annotation = braker_annotation_ch
    miniprot_alignment = miniprot_alignment_ch
    star_pass2 = all_bams_ch
    stringtie_transcripts = stringtie_transcripts_ch
}
