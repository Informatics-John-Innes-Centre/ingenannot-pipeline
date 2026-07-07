
process ingenannot_validate_annotation {
    input:
    tuple val(label), path(annotation)

    output:
    tuple val(label), path(annotation)

    script:
    """
    ingenannot -v 2 validate ${annotation}
    """
}

process sort_bgzip_index_miniprot_mappings {
    input:
    val genome_prefix
    tuple val(label), path(miniprot)

    output:
    tuple path("${genome_prefix}_miniprot_sorted.gff.gz"),
          path("${genome_prefix}_miniprot_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 --parallel=${task.cpus} ${miniprot} > ${genome_prefix}_miniprot_sorted.gff
    bgzip --threads ${task.cpus} ${genome_prefix}_miniprot_sorted.gff
    tabix -C -p gff ${genome_prefix}_miniprot_sorted.gff.gz
    """
}

process combine_stringtie_transcripts {
    input:
    val genome_prefix
    path transcripts

    output:
    path "${genome_prefix}_stringtie_whole.gff"

    script:
    """
    stringtie --merge ${transcripts} -o ${genome_prefix}_stringtie_whole.gff -l ${genome_prefix}_str
    """
}

process sort_bgzip_index_combined_stringtie_transcript {
    input:
    val genome_prefix
    path combined_transcript

    output:
    tuple path("${genome_prefix}_stringtie_whole_sorted.gff.gz"),
          path("${genome_prefix}_stringtie_whole_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 --parallel=${task.cpus} ${combined_transcript} | \
    bgzip -c --threads ${task.cpus} - > ${genome_prefix}_stringtie_whole_sorted.gff.gz
    tabix -C -p gff ${genome_prefix}_stringtie_whole_sorted.gff.gz
    """
}

process create_masked_fasta_minimap_index {
    input:
    val genome_prefix
    path masked_fasta

    output:
    path "${genome_prefix}_masked.fasta.mmi"

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
    val genome_prefix
    path cram_flnc

    output:
    path "${genome_prefix}.flnc.bam"

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
    val genome_prefix
    path masked_fasta_index
    path bam_flnc

    output:
    path "${genome_prefix}_iso.bam"

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
    val genome_prefix
    path aligned_isoseq

    output:
    path "${genome_prefix}_collapse_iso.gff"

    script:
    """
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        ${aligned_isoseq}  \
        ${genome_prefix}_collapse_iso.gff
    """
}

process samtools_index {
    input:
    path target

    output:
    path "${target}.csi"

    script: 
    """
    samtools index -c ${target}
    """
}

process isoform_ranking {

    input:
    val genome_prefix
    path collapsed_isoseq
    path bam_indexes
    path bam_files
    val manifest_str

    output:
    path "${genome_prefix}_isoforms.top.gff"

    script:
    """
    ingenannot -v 2 -p ${task.cpus} isoform_ranking \
        ${collapsed_isoseq} \
        -p ${genome_prefix}_isoforms \
        -f <(printf '%s\\n' "${manifest_str}") \
        --alt_threshold 0.1
    """
}

process sort_bgzip_index_top_isoforms {
    input:
    val genome_prefix
    path top_isoforms

    output:
    tuple path("${genome_prefix}_lr_top_sorted.gff.gz"),
          path("${genome_prefix}_lr_top_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 --parallel=${task.cpus} ${top_isoforms} | \
    bgzip -c --threads ${task.cpus} - > ${genome_prefix}_lr_top_sorted.gff.gz
    tabix -C -p gff ${genome_prefix}_lr_top_sorted.gff.gz
    """
}

process compute_aed_score_for_annotation {
    input:
    val genome_prefix
    tuple val(label), path(annotation)
    tuple path(miniprot_gff), path(miniprot_gff_csi)
    tuple path(stringtie_gff), path(stringtie_gff_csi)
    tuple path(top_isoforms_gff), path(top_isoforms_gff_csi)
    
    output:
    tuple val(label), path("${genome_prefix}_${label}.aed.gff"), path("scatter_hist_aed.${label}.png")

    script:
    """
    ingenannot -v 2 -p ${task.cpus} \
        aed \
        ${annotation} \
        ${genome_prefix}_${label}.aed.gff \
        ${label} \
        ${stringtie_gff} \
        ${miniprot_gff} \
        --longreads ${top_isoforms_gff} \
        --evtrstranded \
        --longreads_source "PacBio" \
        --penalty_overflow 0.2 \
        --aed_tr_cds_only \
        --evtr_source "stringtie" \
        --evpr_source "miniprot"
    """
}

process ingenannot_selection_process {
    input:
    val genome_prefix
    path select_file

    output:
    path("${genome_prefix}_select.genes.gff"), emit: gff
    path("${genome_prefix}_select.genes.gff.scatter_hist_aed.png"), emit: plot

    script:
    """
    ingenannot -v 2 -p ${task.cpus} \
        select \
        ${select_file} \
        ${genome_prefix}_select.genes.gff \
        --noaed \
        --clustranded \
        --nbsrc_filter 2 \
        --aedtr 0.4 \
        --aedpr 0.2 \
        --use_ev_lg \
        --min_cds_len 100 \
        --no_partial \
        --genome ${params.frozDir}/${genome_prefix}/${genome_prefix}.fasta \
        --no_cds_overlap
    """
}

process ingenannot_compare {    
    input:
    path select_output
    path compare_fof

    output:
    path "ingenannot_compare.log"

    script:
    """
    ingenannot -v 2 -p ${task.cpus} \
    compare \
    ${compare_fof} > "ingenannot_compare.log"
    """
}

workflow ingenannot {
    take:
    masked_file
    tiberius_annotation
    helixer_annotation
    braker_annotation
    annevo_annotation
    miniprot_alignment
    star_pass2
    stringtie_transcripts
    genome_prefix
    isoseq_prefix
    
    main:
    
    def cram_flnc = channel.fromPath("${params.isoseqDirectory}/${isoseq_prefix}.flnc.cram")

    // validate the 4 annotations
    def annotations = annevo_annotation
        .mix(braker_annotation)
        .mix(helixer_annotation)
        .mix(tiberius_annotation)
        
    def validated_annotations_ch = ingenannot_validate_annotation(annotations)

    // magic with miniprot
    def miniprot_gff_csi_ch = sort_bgzip_index_miniprot_mappings(genome_prefix, miniprot_alignment)
    
     //combine em to one!
    def combined_stringtie_transcript = combine_stringtie_transcripts(genome_prefix, stringtie_transcripts)
    def stringtie_gff_csi_ch = sort_bgzip_index_combined_stringtie_transcript(genome_prefix, combined_stringtie_transcript)

    def masked_fasta_minimap_index = create_masked_fasta_minimap_index(genome_prefix, masked_file)

    def bam_flnc = convert_cram_flnc_to_bam(genome_prefix, cram_flnc)
    def aligned_isoseq_reads = align_isoseq_reads_to_genome(genome_prefix, masked_fasta_minimap_index, bam_flnc)
    def collapsed_isoseq = collapse_isoseq(genome_prefix, aligned_isoseq_reads)

    ch_manifest_str = star_pass2
        .collect()
        .map { bam_list ->
            bam_list.collect { bam -> "${bam.name}\ttrue\ttrue" }.join('\n')
        }

    def bam_indexes = samtools_index(star_pass2).collect()
    def top_isoforms = isoform_ranking(genome_prefix, collapsed_isoseq, bam_indexes, star_pass2.collect(), ch_manifest_str)   
    def top_isoforms_gff_csi_ch = sort_bgzip_index_top_isoforms(genome_prefix, top_isoforms)

    def aed_scores_ch = compute_aed_score_for_annotation(
        genome_prefix, 
        validated_annotations_ch, 
        miniprot_gff_csi_ch.collect(), 
        stringtie_gff_csi_ch.collect(), 
        top_isoforms_gff_csi_ch.collect()
    )

    
    def select_fof_ch = aed_scores_ch
        .map { annotation_name, aed_gff_path, _aed_scatter_hist ->
            "${aed_gff_path}\t${annotation_name}"
        }
        .collectFile(name: "${genome_prefix}_select.fof", newLine: true, sort: false)

    def select_output = ingenannot_selection_process(genome_prefix, select_fof_ch)


    def compare_fof_ch = select_fof_ch
    .combine(select_output.gff)
    .flatMap { fof, select_path ->
        fof.readLines() + ["${select_path}\tselect"]
    }
    .collectFile(name: "${genome_prefix}_compare.fof", newLine: true, sort: false)

    
    def ingenannot_compare_log_ch = ingenannot_compare(select_output.gff, compare_fof_ch)

    emit: 
    ingenannot_compare_log = ingenannot_compare_log_ch
    aed_scores = aed_scores_ch
    ingenannot_select_gff = select_output.gff
    ingenannot_select_plot = select_output.plot
    miniprot_gff_csi = miniprot_gff_csi_ch
    stringtie_gff_csi = stringtie_gff_csi_ch
    top_isoforms_gff_csi = top_isoforms_gff_csi_ch
}