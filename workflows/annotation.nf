process helixer {
    label 'gpu'
    
    input:
    tuple val(genome_prefix), path(masked_fasta)

    output:
    tuple val(genome_prefix), val("helixer"), path("${genome_prefix}_helixer.gff")

    script:
    """
    mkdir -p tmp
    Helixer.py \
        --downloaded-model-path /opt/helixer/models \
        --lineage ${params.helixerLineage} \
        --fasta-path ${masked_fasta} \
        --species ${genome_prefix} \
        --gff-output-path "${genome_prefix}_helixer.gff" \
        --peak-threshold 0.8 \
        --subsequence-length 64152 \
        --overlap-offset 32076 \
        --overlap-core-length 48114 \
        --temporary-dir tmp
    """
}

process helixer_strip_UTRs {
    input: 
    tuple val(genome_prefix), val("helixer"), path(helixer_gff)

    output:
    tuple val(genome_prefix), val("helixer"), path("${genome_prefix}_helixer_no_UTRs.gff")

    script:
    """
    strip_UTRs.sh --drop-noncoding ${helixer_gff} > "${genome_prefix}_helixer_no_UTRs.gff"
    """
}

process make_tiberius_gene_ids_unique {
    input: 
    tuple val(genome_prefix), val("tiberius"), path(tiberius_gtf)

    output:
    tuple val(genome_prefix), val("tiberius"), path("${genome_prefix}_tiberius_unique_gene_ids.gtf")

    script:
    """
    make_tiberius_gene_ids_unique.sh ${tiberius_gtf} > ${genome_prefix}_tiberius_unique_gene_ids.gtf
    """
}

process convert_tiberius_gtf_to_gff {
    input:
    tuple val(genome_prefix), val("tiberius"), path(tiberius_gtf)

    output:
    tuple val(genome_prefix), val("tiberius"), path("${genome_prefix}_tiberius.gff")

    script:
    """
    agat_convert_sp_gxf2gxf.pl \
        -gff ${tiberius_gtf} \
        -o "${genome_prefix}_tiberius.gff" \
        -gvo 3 \
        --no_progressbar
    """
}

process tiberius {
    label 'gpu'

    input:
    tuple val(genome_prefix), path(masked_fasta)

    output:
    tuple val(genome_prefix), val("tiberius"), path("${genome_prefix}_tiberius.gtf")

    script:
    """
    python3 /opt/Tiberius/tiberius.py \
        --genome ${masked_fasta} \
        --model /opt/Tiberius/model_weights/${params.tiberiusModel} \
        --out ${genome_prefix}_tiberius.gtf
    """
}

process annevo {
    label 'gpu'

    input:
    tuple val(genome_prefix), path(masked_fasta)

    output:
    tuple val(genome_prefix), val("annevo"), path("${genome_prefix}_annevo.gff")

    script:
    """
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"
    mkdir -p "\$NUMBA_CACHE_DIR"

    annotation.py \
        --genome ${masked_fasta} \
        --model_path /usr/local/share/annevo-2.3.2/saved_model/${params.annevoModel} \
        --output ${genome_prefix}_annevo.gff \
        --lineage ${params.annevoLineage} \
        --threads ${task.cpus} \
        --num_workers ${task.cpus}
    """
}

process braker3 {    
    input:
    tuple val(genome_prefix), path(masked_fasta), path(proteinDatabase), path(bam_files)

    output:
    tuple val(genome_prefix), val("braker3"), path("${genome_prefix}_braker3.gff")

    script:
    """
    cp -r /opt/Augustus/config ./augustus_config
    export AUGUSTUS_CONFIG_PATH=\$PWD/augustus_config
    braker.pl \
    --prot_seq=${proteinDatabase} \
    --genome=${masked_fasta} \
    --bam=${bam_files.join(',')} \
    --threads=${task.cpus} \
    --species=${genome_prefix} \
    --workingdir=\$PWD \
    --gff3 \
    --verbosity=4
    mv braker.gff3 ${genome_prefix}_braker3.gff
    """
}


process ingenannot_validate_annotation {
    input:
    tuple val(genome_prefix), val(label), path(annotation)

    output:
    tuple val(genome_prefix), val(label), path(annotation)

    script:
    """
    ingenannot -v 2 validate ${annotation}
    """
}

workflow annotate {
    take:
    masked_file
    protein_databases
    bam_files

    main:
    def helixer_annotation_ch = helixer(masked_file)
    def helixer_no_UTRs = helixer_strip_UTRs(helixer_annotation_ch)

    def tiberius_annotation_ch = tiberius(masked_file)
    def tiberius_unique_gene_ids = make_tiberius_gene_ids_unique(tiberius_annotation_ch)
    def tiberius_gff = convert_tiberius_gtf_to_gff(tiberius_unique_gene_ids)
    def annevo_annotation_ch = annevo(masked_file)

    def braker3_annotation_ch = braker3(masked_file.join(protein_databases).join(bam_files))

    def all_annotations = helixer_no_UTRs
        .mix(tiberius_gff)
        .mix(annevo_annotation_ch)
        .mix(braker3_annotation_ch)

    def validated_annotations_ch = ingenannot_validate_annotation(all_annotations)

    emit:
    annotations = validated_annotations_ch
}