process create_masked_fasta_minimap_index {
    input:
    tuple val(genome_prefix), path(masked_fasta)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_masked.fasta.mmi")

    script:
    """
    pbmm2 index \
        --preset ISOSEQ \
        ${masked_fasta} \
        ${genome_prefix}_masked.fasta.mmi
    """
}

process convert_cram_flnc_to_bam {
    input:
    tuple val(genome_prefix), path(cram_flnc)

    output:
    tuple val(genome_prefix), path("${genome_prefix}.flnc.bam")

    script:
    """
    samtools view \
    -b  \
    -o ${genome_prefix}.flnc.bam \
    ${cram_flnc}
    """
}

process align_isoseq_reads_to_genome {
    input:
    tuple val(genome_prefix), path(bam_flnc), path(masked_fasta_index)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_iso.bam")

    script:
    """
    pbmm2 align \
        ${masked_fasta_index} \
        ${bam_flnc}\
        ${genome_prefix}_iso.bam \
        --sort \
        -j ${task.cpus} \
        -J ${task.cpus} \
        --bam-index CSI
    """
}

process collapse_isoseq {
    input: 
    tuple val(genome_prefix), path(aligned_isoseq)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_collapse_iso.gff")

    script:
    """
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        ${aligned_isoseq}  \
        ${genome_prefix}_collapse_iso.gff
    """
}

workflow isoseq {
    take: 
    masked_file
    cram_flnc

    main:
    def masked_fasta_minimap_index = create_masked_fasta_minimap_index(masked_file)

    def bam_flnc = convert_cram_flnc_to_bam(cram_flnc)

    def joined = bam_flnc.join(masked_fasta_minimap_index)

    def aligned_isoseq = align_isoseq_reads_to_genome(joined)
    
    def collapsed_isoseq_ch = collapse_isoseq(aligned_isoseq)

    emit:
    collapsed_isoseq = collapsed_isoseq_ch
}