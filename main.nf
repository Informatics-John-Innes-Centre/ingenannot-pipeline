#!/usr/bin/env nextflow

process softmask {    
    input:
    val genome_prefix

    output:
    path "${genome_prefix}_masked.fasta"

    script:
    """
    bedtools maskfasta \
        -soft \
        -fi ${params.frozDir}/${genome_prefix}/${genome_prefix}.fasta \
        -bed ${params.frozDir}/${genome_prefix}/${genome_prefix}_all_repeats.bed \
        -fo ${genome_prefix}_masked.fasta
    """
}

params {
    frozDir
    helixerModelDir
    annevoModelDir
    tiberiusModelDir
    accessionFile
    proteinDatabase
    fastqDirectory
    isoseqDirectory
    lineNumber 
}

include { annotate } from './workflows/annotate.nf'
include { ingenannot } from './workflows/ingenannot.nf'

workflow {
    main:
    def accessions = file(params.accessionFile).readLines()
    def target = accessions[params.lineNumber as int]
    def accession = target.split('\t')[0]
    def genome_prefix = target.split('\t')[1]
    def illumina_prefix = target.split('\t')[2]
    def isoseq_prefix = target.split('\t')[3]

    println "Accession: ${accession}"
    println "Prefix: ${genome_prefix}"
    println "Illumina RNASEQ Prefix: ${illumina_prefix}"
    println "ISOSEQ Prefix: ${isoseq_prefix}"

    def masked_file = softmask(genome_prefix)

    def annotate_result = annotate(masked_file, genome_prefix, illumina_prefix)
    def ingenannot_result = ingenannot(
        masked_file, 
        annotate_result.tiberius_annotation, 
        annotate_result.helixer_annotation, 
        annotate_result.braker_annotation, 
        annotate_result.annevo_annotation, 
        annotate_result.miniprot_alignment, 
        annotate_result.star_pass2, 
        annotate_result.stringtie_transcripts, 
        genome_prefix, 
        isoseq_prefix
    )
    //  HARCODED
    /* 
    def tiberius_annotation = channel.of(tuple('tiberius', file("/jic/scratch/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/JI1006_v3.0_tiberius_2.0.6.gtf")))    
    def miniprot_alignment = channel.of(tuple('miniprot', file("/jic/scratch/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/JI1006_v3.0_miniprot.gff3")))
    def annevo_annotation = channel.of(tuple('annevo', file("/jic/scratch
/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/JI1006_v3.0_annevo.gff")))
    def braker_annotation = channel.of(tuple('braker', file("/jic/scratch/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/JI1006_v3.0_braker3.gff")))
    def helixer_annotation = channel.of(tuple('helixer', file("/jic/scratch/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/JI1006_v3.0_helixer.gff")))
    def all_bams = channel.fromPath("/jic/scratch/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/star_pass2/*")
    channel
        .fromPath("/jic/scratch/platforms/informatics/peaterpan/nextflow_annotation/inputs/ingenannot_inputs/stringtie/Piful_JI1006-*.gtf")
        .collect()
        .set { stringtie_transcripts }
    def ingenannot_result = ingenannot(masked_file, tiberius_annotation, helixer_annotation, braker_annotation, annevo_annotation, miniprot_alignment, all_bams, stringtie_transcripts, genome_prefix, isoseq_prefix, )
    */
    // HARDCODED
    publish:
    tiberius_annotation = annotate_result.tiberius_annotation
    annevo_annotation = annotate_result.annevo_annotation
    braker_annotation = annotate_result.braker_annotation
    helixer_annotation = annotate_result.helixer_annotation
    ingenannot_compare_log = ingenannot_result.ingenannot_compare_log
    aed_scores = ingenannot_result.aed_scores
    star_pass2 = annotate_result.star_pass2
    stringtie_transcripts = annotate_result.stringtie_transcripts
    ingenannot_select_gff  = ingenannot_result.ingenannot_select_gff
    ingenannot_select_plot = ingenannot_result.ingenannot_select_plot
    miniprot_gff_csi = ingenannot_result.miniprot_gff_csi
    stringtie_gff_csi = ingenannot_result.stringtie_gff_csi
    top_isoforms_gff_csi = ingenannot_result.top_isoforms_gff_csi
}

output {
    ingenannot_compare_log {}
    // aed scores
    aed_scores {
        path { file -> "aed_scores/${file.name}" }
    }
    // annotations
    tiberius_annotation {
        path { _id, file -> "annotations/${file.name}" }
    }
    annevo_annotation {
        path { _id, file -> "annotations/${file.name}" }
    }
    braker_annotation {
        path { _id, file -> "annotations/${file.name}" }
    }
    helixer_annotation {
        path { _id, file -> "annotations/${file.name}" }
    }
    star_pass2 {
        path { _id, file -> "star_pass2/${file.name}"}
    }
    stringtie_transcripts {
        path { file -> "stringtie/${file.name}"}
    }
    ingenannot_select_gff {}
    ingenannot_select_plot {}
    miniprot_gff_csi {}
    stringtie_gff_csi {}
    top_isoforms_gff_csi {}
}
