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
        --lineage land_plant \
        --fasta-path ${masked_fasta} \
        --species Pisum \
        --gff-output-path "${genome_prefix}_helixer.gff" \
        --peak-threshold 0.8 \
        --subsequence-length 64152 \
        --overlap-offset 32076 \
        --overlap-core-length 48114 \
        --temporary-dir tmp
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
    python /opt/Tiberius/tiberius.py \
        --genome ${masked_fasta} \
        --model /opt/Tiberius/model_weights/angiosperms_softmasking \
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

    python /opt/ANNEVO/annotation.py \
        --genome ${masked_fasta} \
        --model_path /opt/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt \
        --output ${genome_prefix}_annevo.gff \
        --lineage Magnoliopsida \
        --threads ${task.cpus} \
        --num_workers ${task.cpus}
    """
}

process braker3 {    
    input:
    tuple val(genome_prefix), path(masked_fasta), path(bam_files)

    output:
    tuple val(genome_prefix), val("braker3"), path("braker.gff3")

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
    bam_files

    main:
    def helixer_annotation_ch = helixer(masked_file)
    def tiberius_annotation_ch = tiberius(masked_file)
    def annevo_annotation_ch = annevo(masked_file)

    def braker3_annotation_ch = braker3(masked_file.join(bam_files))

    def all_annotations = helixer_annotation_ch
        .mix(tiberius_annotation_ch)
        .mix(annevo_annotation_ch)
        .mix(braker3_annotation_ch)

    def validated_annotations_ch = ingenannot_validate_annotation(all_annotations)

    emit:
    annotations = validated_annotations_ch
}