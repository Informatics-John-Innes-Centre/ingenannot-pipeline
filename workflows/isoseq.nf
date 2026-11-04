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

process create_pbi_index_for_bam {
    input:
    tuple val(genome_prefix), path(bam)

    output:
    tuple val(genome_prefix), path("${bam}.pbi")

    script: 
    """
    pbindex ${bam}
    """
}

process convert_flnc_bam_to_clustered_bam {
    input:
    tuple val(genome_prefix), path(bam_flnc)

    output:
    tuple val(genome_prefix), path("${genome_prefix}.clustered.bam"), path("${genome_prefix}.clustered.bam.pbi")

    script:
    """
    isoseq cluster2 \
        -j ${task.cpus} \
        ${bam_flnc} \
        ${genome_prefix}.clustered.bam
    """
}

process align_isoseq_reads_to_genome {
    input:
    tuple val(genome_prefix), path(bam), path(masked_fasta_index)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_clustered_aligned.bam")

    script:
    """
    pbmm2 align \
        ${masked_fasta_index} \
        ${bam}\
        ${genome_prefix}_clustered_aligned.bam \
        --sort \
        -j ${task.cpus} \
        -J ${task.cpus} \
        --bam-index CSI
    """
}

process collapse_isoseq {
    input: 
    tuple val(genome_prefix), path(clustered_bam), path(clustered_bam_index),  path(flnc_bam), path(flnc_bam_pbi)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_collapse_iso.gff"), path("${genome_prefix}_collapse_iso.abundance.txt")

    script:
    """
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        ${clustered_bam}  \
        ${flnc_bam} \
        ${genome_prefix}_collapse_iso.gff
    """
}

process sort_bgzip_index_collapsed_isoseq {
    input:
    tuple val(genome_prefix), path(collapsed_isoseq)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_lr_raw_sorted.gff.gz"), path("${genome_prefix}_lr_raw_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 \
        ${collapsed_isoseq} | \
        bgzip -c --threads ${task.cpus} - > ${genome_prefix}_lr_raw_sorted.gff.gz
    tabix -C -p gff ${genome_prefix}_lr_raw_sorted.gff.gz
    """
}

workflow isoseq {
    take: 
    masked_file
    cram_flnc

    main:
    def masked_fasta_minimap_index = create_masked_fasta_minimap_index(masked_file)

    def bam_flnc = convert_cram_flnc_to_bam(cram_flnc)
    def bam_flnc_pbi_index = create_pbi_index_for_bam(bam_flnc)
    def clustered_bam_and_index = convert_flnc_bam_to_clustered_bam(bam_flnc)
    
    // alignment doesn't need the pbi index
    def clustered_bam = clustered_bam_and_index.map { genome_prefix, bam, _index -> tuple(genome_prefix, bam )}
    def joined = clustered_bam.join(masked_fasta_minimap_index)
    def aligned_isoseq = align_isoseq_reads_to_genome(joined)
    
    // collapse needs a pbi index of the bam_flnc and the pbi index from the clustered bam
    def clustered_bam_index = clustered_bam_and_index.map {genome_prefix, _bam, index -> tuple(genome_prefix, index)}
    def collapse_input = aligned_isoseq.join(clustered_bam_index).join(bam_flnc).join(bam_flnc_pbi_index)
    def collapse_isoseq_result = collapse_isoseq(collapse_input)
    def collapsed_isoseq_ch = collapse_isoseq_result.map { genome_prefix, collapse_isoseq, _abundance -> tuple(genome_prefix, collapse_isoseq)}
    
    def sorted_isoseq_ch = sort_bgzip_index_collapsed_isoseq(collapsed_isoseq_ch)
    def abundance_txt_ch = collapse_isoseq_result.map { genome_prefix, _collapse_isoseq, abundance -> tuple(genome_prefix, abundance)}

    emit:
    sorted_isoseq = sorted_isoseq_ch
    collapsed_isoseq = collapsed_isoseq_ch
    abundance_txt = abundance_txt_ch
    long_read_bams = aligned_isoseq
}