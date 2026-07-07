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
