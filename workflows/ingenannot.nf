process isoform_ranking {
    input:
    tuple val(genome_prefix), path(collapsed_isoseq), path(bams), path(indices), val(manifest_str)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_isoforms.top.gff")

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
    tuple val(genome_prefix), path(top_isoforms)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_lr_top_sorted.gff.gz"), path("${genome_prefix}_lr_top_sorted.gff.gz.csi")

    script:
    """
    sort -k1,1 -k4g,4 --parallel=${task.cpus} ${top_isoforms} | \
    bgzip -c --threads ${task.cpus} - > ${genome_prefix}_lr_top_sorted.gff.gz
    tabix -C -p gff ${genome_prefix}_lr_top_sorted.gff.gz
    """
}


process compute_aed_score_for_annotation {
    input:
    tuple val(genome_prefix), val(label), path(annotation), path(miniprot_gff), path(miniprot_gff_csi), path(stringtie_gff), path(stringtie_gff_csi), path(top_isoforms_gff), path(top_isoforms_gff_csi)

    output:
    tuple val(genome_prefix), val(label), path("${genome_prefix}_${label}.aed.gff"), path("scatter_hist_aed.${label}.png")

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
    tuple val(genome_prefix), val(select_file_contents)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_select.genes.gff"), path("${genome_prefix}_select.genes.gff.scatter_hist_aed.png")
    script:
    """
    printf '%s\n' "${select_file_contents}" > ${genome_prefix}_select.fof
    ingenannot -v 2 -p ${task.cpus} \
        select \
        ${genome_prefix}_select.fof \
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
    tuple val(genome_prefix), val(compare_file_contents)

    output:
    tuple val(genome_prefix), path("ingenannot_compare.log")

    script:
    """
    printf '%s\n' "${compare_file_contents}" > ${genome_prefix}_compare.fof
    ingenannot -v 2 -p ${task.cpus} \
    compare \
    ${genome_prefix}_compare.fof > "ingenannot_compare.log"
    """
}

workflow ingenannot {
    take:
    annotations
    collapsed_isoseq
    bam_files
    bam_indices
    miniprot_gff_csi
    stringtie_gff_csi

    main:

    def manifest = bam_files
        .map { genome, bams ->
            tuple(
                genome,
                bams.collect { "${it.name}\ttrue\ttrue" }.join('\n')
            )
        }
    def isoform_ranking_input = collapsed_isoseq
        .join(bam_files)
        .join(bam_indices)
        .join(manifest)

    def top_isoforms = isoform_ranking(isoform_ranking_input)

    def top_isoforms_gff_csi_ch =
        sort_bgzip_index_top_isoforms(top_isoforms)

    def aed_input = annotations
        .combine(miniprot_gff_csi, by: 0)
        .combine(stringtie_gff_csi, by: 0)
        .combine(top_isoforms_gff_csi_ch, by: 0)

    def aed_scores_ch =
        compute_aed_score_for_annotation(aed_input)

    def select_fof_input = aed_scores_ch
        .map { genome_prefix, label, aed_gff, _plot ->
            tuple(groupKey(genome_prefix, 4), "${aed_gff}\t${label}")
        }
        .groupTuple()
        .map { gkey, lines -> tuple(gkey.getGroupTarget(), lines) }

    def selection_process_input = select_fof_input
        .map { genome_prefix, lines ->
            tuple(genome_prefix, lines.join('\n'))
        }

    def select_output =
        ingenannot_selection_process(selection_process_input)

    def comparison_process_input = select_fof_input
        .map { genome_prefix, lines ->
            tuple(
                genome_prefix,
                (lines + "${genome_prefix}_select.genes.gff\tselect")
                    .join('\n')
            )
        }

    def ingenannot_compare_log_ch =
        ingenannot_compare(comparison_process_input)

    emit:
    top_isoforms_gff_csi = top_isoforms_gff_csi_ch
    aed_scores = aed_scores_ch
    selection = select_output
    comparison_log = ingenannot_compare_log_ch
}

