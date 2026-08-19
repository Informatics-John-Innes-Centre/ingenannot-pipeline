process miniprot_alignment  {
    input:
    tuple val(genome_prefix), path(masked_fasta), path(proteinDatabase)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_miniprot.gff3")

    script:
    """
    miniprot \
        -t ${task.cpus} \
        -d ${genome_prefix}_masked.mpi \
        ${masked_fasta}

    miniprot \
        -t ${task.cpus} \
        ${genome_prefix}_masked.mpi  \
        ${proteinDatabase} \
        --gff \
        > ${genome_prefix}_miniprot.gff3
    """
}

process sort_bgzip_index_miniprot_mappings {
    input:
    tuple val(genome_prefix), path(miniprot)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_miniprot_sorted.gff.gz"), path("${genome_prefix}_miniprot_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 --parallel=${task.cpus} ${miniprot} > ${genome_prefix}_miniprot_sorted.gff
    bgzip --threads ${task.cpus} ${genome_prefix}_miniprot_sorted.gff
    tabix -C -p gff ${genome_prefix}_miniprot_sorted.gff.gz
    """
}

workflow miniprot {
    take:
    masked_fasta
    protein_databases

    main:
    def miniprot_input = masked_fasta.join(protein_databases)
    def miniprot_ch = miniprot_alignment(miniprot_input)
    def sorted_miniprot_gff_csi = sort_bgzip_index_miniprot_mappings(miniprot_ch)
    
    emit:
    gff_csi = sorted_miniprot_gff_csi
}