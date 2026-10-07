/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_AFTER } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'
include { GFFREAD             } from '../modules/nf-core/gffread/main'
include { STAR_ALIGN   } from '../modules/nf-core/star/align/main'
include { SALMON_QUANT } from '../modules/nf-core/salmon/quant/main'
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'
include { CUSTOM_TX2GENE } from '../modules/nf-core/custom/tx2gene/main'
include { TXIMETA_TXIMPORT } from '../modules/nf-core/tximeta/tximport/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir
    fasta
    gtf

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    //
    // MODULES: Prepare reference index and transcript sequences
    //
    if (!fasta || !gtf) {
        error 'Provide both --fasta and --gtf'
    }

    def fasta_file = file(fasta, checkIfExists: true)
    def gtf_file   = file(gtf, checkIfExists: true)

    def ch_fasta = channel.value(
        [[id: 'reference'], fasta_file]
    )

    def ch_gtf = channel.value(
        [[id: 'reference'], gtf_file]
    )

    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
    //
    // Prepare the transcript reference for Salmon
    //
    def ch_salmon_reference

    if (params.transcript_fasta) {

        def transcript_fasta_file = file(
            params.transcript_fasta,
            checkIfExists: true
        )

        ch_salmon_reference = channel.value([
            [id: 'reference'],
            [],
            gtf_file,
            transcript_fasta_file
        ])

    } else {

        GFFREAD(
            ch_gtf,
            channel.value(fasta_file)
        )

        ch_salmon_reference = GFFREAD.out.gffread_fasta
            .map { meta, transcript_fasta ->
                [meta, [], gtf_file, transcript_fasta]
            }
            .first()
    }
    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'rnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    //
    // MODULES: Select trimming tool and run post-trimming QC
    //
    def ch_trimmed_reads

    if (params.trimmer == 'trimgalore') {

        TRIMGALORE(ch_samplesheet)

        ch_trimmed_reads = TRIMGALORE.out.reads

        ch_multiqc_files = ch_multiqc_files.mix(
            TRIMGALORE.out.log.map { meta, logs -> logs }
        )

    } else if (params.trimmer == 'fastp') {

        // FASTP requires [meta, reads, adapter_fasta].
        // An empty list means no custom adapter FASTA is provided.
        def ch_fastp_input = ch_samplesheet.map { meta, reads ->
            [meta, reads, []]
        }

        FASTP(
            ch_fastp_input,
            false,  // Keep passing trimmed reads
            false,  // Do not save failed reads
            false   // Do not merge paired reads
        )

        ch_trimmed_reads = FASTP.out.reads

        ch_multiqc_files = ch_multiqc_files.mix(
            FASTP.out.json.map { meta, report -> report }
        )

    } else {
        error "Invalid --trimmer '${params.trimmer}'. Choose 'trimgalore' or 'fastp'."
    }

    FASTQC_AFTER(ch_trimmed_reads)

    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_AFTER.out.zip.map { meta, files -> files }
    )
    //
    // MODULE: Align trimmed reads with STAR
    //
    def ch_star_index = STAR_GENOMEGENERATE.out.index.first() // for multiple samples, only take the first index (they are all the same)
    
    STAR_ALIGN(
        ch_trimmed_reads,
        ch_star_index,
        ch_gtf,
        true
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        STAR_ALIGN.out.log_final.map { meta, log -> log }
    )


    //
    // MODULES: Sort genome alignments and mark duplicates
    //
    def ch_bam_reference = channel.value(
        [[id: 'reference'], [], []]
    )

    PICARD_MARKDUPLICATES(
        STAR_ALIGN.out.bam_sorted_aligned,
        ch_bam_reference
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        PICARD_MARKDUPLICATES.out.metrics.map { meta, metrics -> metrics }
    )
    //
    // MODULE: Quantify transcript alignments with Salmon
    //
    SALMON_QUANT(
        STAR_ALIGN.out.bam_transcript,
        ch_salmon_reference
    )

    ch_multiqc_files = ch_multiqc_files.mix(
        SALMON_QUANT.out.results.map { meta, directory -> directory }
    )
    //
    // MODULE: Build transcript-to-gene mapping
    //
    def ch_tx2gene_quant = SALMON_QUANT.out.results
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .filter { samples -> !samples.isEmpty() }
        .map { samples -> samples.first() }

    CUSTOM_TX2GENE(
        ch_gtf,
        ch_tx2gene_quant,
        'salmon',
        'gene_id',
        'gene_name'
    )
    //
    // MODULE: Merge quantifications across samples
    //
    def ch_all_quants = SALMON_QUANT.out.results
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .filter { samples -> !samples.isEmpty() }
        .map { samples ->
            [
                [id: 'salmon_merged'],
                samples.collect { sample -> sample[1] }
            ]
        }

    TXIMETA_TXIMPORT(
        ch_all_quants,
        CUSTOM_TX2GENE.out.tx2gene,
        'salmon'
    )
    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'rnaseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    gene_tpm = TXIMETA_TXIMPORT.out.tpm_gene
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
