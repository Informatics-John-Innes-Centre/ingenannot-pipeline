process stringtie_assemble {
    input:
    tuple val(genome_prefix), val(sample_id), path(bam)

    output:
    tuple val(genome_prefix), path("${sample_id}.gtf")

    script:"""
    stringtie \
        -p ${task.cpus} \
        --nasc \
        ${bam} \
        > ${sample_id}.gtf \

    """
}

process combine_stringtie_transcripts {
    input:
    tuple val(genome_prefix), path(transcripts)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_stringtie_whole.gff")

    script:
    """
    stringtie --merge ${transcripts} -o ${genome_prefix}_stringtie_whole.gff -l ${genome_prefix}_str
    """
}

process sort_bgzip_index_combined_stringtie_transcript {
    input:
    tuple val(genome_prefix), path(combined_transcript)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_stringtie_whole_sorted.gff.gz"), path("${genome_prefix}_stringtie_whole_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 --parallel=${task.cpus} ${combined_transcript} | \
    bgzip -c --threads ${task.cpus} - > ${genome_prefix}_stringtie_whole_sorted.gff.gz
    tabix -C -p gff ${genome_prefix}_stringtie_whole_sorted.gff.gz
    """
}

workflow stringtie {
    take: 
    star_bams

    main:
    def bams = star_bams.map { genome, sample, bam, _csi ->
        tuple(genome, sample, bam)
    }

    def transcripts = stringtie_assemble(bams)
    def grouped_transcripts = transcripts.groupTuple()
    def merged = combine_stringtie_transcripts(grouped_transcripts) 
    def sorted = sort_bgzip_index_combined_stringtie_transcript(merged)

    emit:
    gff_csi = sorted
}