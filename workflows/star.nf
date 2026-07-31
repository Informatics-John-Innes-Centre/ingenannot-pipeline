process star_index {
    input:
    tuple val(genome_prefix), path(masked_fasta)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_star_index/")

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
    tuple val(genome_prefix), val(sample_id), path(reads)

    output:
    tuple val(genome_prefix), val(sample_id), path(["trimmed-${reads[0]}", "trimmed-${reads[1]}"])

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
    tuple val(genome_prefix), path(star_index), val(sample_id), path(reads)

    output:
    tuple val(genome_prefix), val(sample_id), path("${sample_id}_SJ.out.tab")

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
    tuple val(genome_prefix), path(star_index), val(sample_id), path(reads), path(sj_files)

    output:
    tuple val(genome_prefix), val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

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

process samtools_index {
    input:
    tuple val(genome_prefix), val(sample), path(target)

    output:
    tuple val(genome_prefix), val(sample), path(target), path("${target}.csi")

    script: 
    """
    samtools index -c ${target}
    """
}

workflow star {
    take:
    masked_fasta
    rnaseq_pairs

    main:
    def star_index_ch = star_index(masked_fasta)

    def trimmed_rnaseq_ch = trim_rnaseq_data(rnaseq_pairs)

    def first_star_pass_input = star_index_ch.combine(trimmed_rnaseq_ch, by: 0)

    def first_pass_mapped = star_map_rnaseq_data_first_pass(first_star_pass_input)

    // Collect SJ files belonging to each genome
    def all_sj = first_pass_mapped
        .map { genome_prefix, sample_id, sj_file ->
            tuple(genome_prefix, sj_file)
        }
        .groupTuple()
        .map { genome_prefix, sj_files ->
            tuple(genome_prefix, sj_files.sort { it.name })
        }

    def second_pass_input = first_star_pass_input.combine(all_sj, by: 0)

    def second_star_pass = star_map_rnaseq_data_second_pass(second_pass_input)

    // Index every BAM while preserving metadata
    def indexed_bams = samtools_index(second_star_pass)

    emit:
    bams = indexed_bams
}