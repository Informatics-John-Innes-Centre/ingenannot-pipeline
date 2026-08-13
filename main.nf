#!/usr/bin/env nextflow

include { annotate } from './workflows/annotation.nf'
include { miniprot } from './workflows/miniprot.nf'
include { stringtie } from './workflows/stringtie.nf'
include { isoseq } from './workflows/isoseq.nf'
include { star } from './workflows/star.nf'
include { ingenannot } from './workflows/ingenannot.nf'



process softmask {
    input:
    tuple val(genome_prefix), path(genome), path(repeats)

    output:
    tuple val(genome_prefix), path("${genome_prefix}_masked.fasta")

    script:
    """
    bedtools maskfasta \
        -soft \
        -fi ${genome} \
        -bed ${repeats} \
        -fo ${genome_prefix}_masked.fasta
    """
}

params {
    // data inputs
    frozDir
    accessionFile
    proteinDatabase
    fastqDirectory
    isoseqDirectory

    // tool specific configurations
    helixerLineage
    tiberiusModel
    annevoModel
    annevoLineage
}

workflow genome {
    take:
    genome_prefixes
    cram_flnc
    rnaseq_pairs
    sample_counts

    main:
    def softmasking_input = genome_prefixes.map{ genome_prefix -> 
        tuple(
            genome_prefix, 
            file("${params.frozDir}/${genome_prefix}/${genome_prefix}.fasta"), 
            file("${params.frozDir}/${genome_prefix}/${genome_prefix}_all_repeats.bed"))
    }
    def masked_file = softmask(softmasking_input)    
    def star = star(masked_file, rnaseq_pairs, sample_counts)

    def bam_files = star.bams
        .map { genome, sample, bam, _csi ->
            tuple(genome, tuple(sample, bam))
        }
        .combine(sample_counts, by: 0)
        .map { genome, entry, count ->
            tuple(groupKey(genome, count), entry)
        }
        .groupTuple()
        .map { gkey, entries ->
            tuple(
                gkey.toString(),
                entries.sort { entry -> entry[0] }.collect { entry -> entry[1] }
            )
        }

    def bam_indices = star.bams
        .map { genome, sample, _bam, csi ->
            tuple(genome, tuple(sample, csi))
        }
        .combine(sample_counts, by: 0)
        .map { genome, entry, count ->
            tuple(groupKey(genome, count), entry)
        }
        .groupTuple()
        .map { gkey, entries ->
            tuple(
                gkey.toString(),
                entries.sort { entry -> entry[0] }.collect { entry -> entry[1] }
            )
        }

    def miniprot_result = miniprot(masked_file)
    def annotate_result = annotate(masked_file, bam_files)
    def stringtie_result = stringtie(star.bams, sample_counts)
    def isoseq_result = isoseq(masked_file, cram_flnc)
    def ingenannot_result = ingenannot(annotate_result.annotations, isoseq_result.collapsed_isoseq, bam_files, bam_indices, miniprot_result.gff_csi, stringtie_result.gff_csi)

    
    emit:
    masked_genomes      = masked_file
    annotations = annotate_result.annotations
    star_bams = star.bams
    miniprot = miniprot_result.gff_csi
    stringtie = stringtie_result.gff_csi
    aed_scores = ingenannot_result.aed_scores
    top_isoforms_gff_csi = ingenannot_result.top_isoforms_gff_csi
    ingenannot_selection = ingenannot_result.selection
    ingenannot_comparison = ingenannot_result.comparison_log
}

workflow {
    main:
    // READ the csv file
    channel
        .fromPath(params.accessionFile)
        .splitCsv(header: true)
        .map { r ->
            tuple(
                r.genome_prefix,
                r.illumina_prefix,
                r.isoseq_prefix
            )
        }
        .set { rows }

    def genome_prefixes = rows.map { genome_prefix, _illumina_prefix, _isoseq_prefix -> genome_prefix }

    // get isoseq cram files
    def cram_flnc = rows
        .map { genome_prefix, _illumina_prefix, isoseq_prefix ->
            tuple(
                genome_prefix,
                file("${params.isoseqDirectory}/${isoseq_prefix}.flnc.cram")
            )
        }

    // Find the rnaseq pairs
    def rnaseq_files = channel
        .fromFilePairs("${params.fastqDirectory}/*-r{1,2}.fastq.gz")
        .map { sample_id, reads ->
            tuple(sample_id, reads)
        }
    
    // Filter out the rnaseq pairs for our illumina prefix
    def rnaseq_pairs = rows
        .map { genome_prefix, illumina_prefix, _isoseq_prefix ->
            tuple(genome_prefix, illumina_prefix)
        }
        .combine(rnaseq_files)
        .filter { _genome_prefix, illumina_prefix, sample_id, _reads ->
            sample_id.startsWith(illumina_prefix)
        }
        .map { genome_prefix, _illumina_prefix, sample_id, reads ->
            tuple(genome_prefix, sample_id, reads)
        }

    def sample_counts = rnaseq_pairs.groupTuple().map { genome_prefix, _sample_id, reads ->
        tuple(genome_prefix, reads.size())
    }

    sample_counts.view()

    def genome_result_ch = genome(genome_prefixes, cram_flnc, rnaseq_pairs, sample_counts)

    publish:
    masked_genomes = genome_result_ch.masked_genomes
    annotations = genome_result_ch.annotations
    star_bams = genome_result_ch.star_bams
    miniprot = genome_result_ch.miniprot
    stringtie = genome_result_ch.stringtie
    aed_scores = genome_result_ch.aed_scores
    top_isoforms_gff_csi = genome_result_ch.top_isoforms_gff_csi
    ingenannot_selection = genome_result_ch.ingenannot_selection
    ingenannot_comparison = genome_result_ch.ingenannot_comparison
}

output {
    masked_genomes {
        path { genome_prefix, _file ->
            "${genome_prefix}"
        }
    }

    annotations {
        path { genome_prefix, _label, _file ->
            "${genome_prefix}/annotations"
        }
    }

    star_bams {
        path { genome_prefix, _sample, _bam, _csi ->
            "${genome_prefix}/star_mappings"
        }
    }

    miniprot {
        path { genome_prefix, _gff, _csi ->
            "${genome_prefix}/miniprot"            
        }
    }

    stringtie {
        path { genome_prefix, _gff, _csi ->
            "${genome_prefix}/stringtie"
        }
    }

    top_isoforms_gff_csi {
        path { genome_prefix, _gff, _csi ->
            "${genome_prefix}/top_isoforms"
        }
    }

    aed_scores {
        path { genome_prefix, _label, _aed_gff, _plot ->
            "${genome_prefix}/ingenannot"
        }
    }

    ingenannot_selection {
        path { genome_prefix, _select_gff, _histogram ->
            "${genome_prefix}/ingenannot"
        }
    }

    ingenannot_comparison {
        path { genome_prefix, _log ->
            "${genome_prefix}/ingenannot"
        }
    }
}