/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowZerocaller.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.intervals ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

if (!["all", "variant-calling", "annotation"].contains(params.stage.toLowerCase())){

    exit 1, "Stage parameter ('${params.stage}') must be one of the following values ['all', 'variant-calling', 'annotation']!"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { VARIANT_CALLING } from '../subworkflows/local/variant_calling.nf'
include {ANNOTATE} from '../subworkflows/local/annotate.nf'

include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'


workflow ZEROCALLER {

    ch_versions = Channel.empty()

    if (["all", "variant-calling"].contains(params.stage.toLowerCase())){
        VARIANT_CALLING(
            params.input,
            params.fasta,
            params.intervals,
            params.fastq_max_size
        )

        ch_versions.mix(VARIANT_CALLING.out.versions).set{ch_versions}

        ch_in_annotate_vcf = VARIANT_CALLING.out.gvcf_parts

    }else{
        Channel
            .fromPath(params.input, checkIfExists:true )
            .splitCsv(header: true)
            .map {
                row -> {
                if (row.tbi)
                    {
                        [row.sample_id,
                        file(row.vcf, checkIfExists:true), 
                        file(row.tbi, checkIfExists:true)]
                    }else{
                        [row.sample_id, file(row.vcf, checkIfExists:true)]
                    }
                }
            }
            .set{ch_in_annotate_vcf}

    }


    if (["all", "annotation"].contains(params.stage.toLowerCase())){
        ANNOTATE(
            ch_in_annotate_vcf,
            params.fasta,
            params.vep_species,
            params.vep_cache_version,
            params.vep_cache,
            params.vep_plugins,
            params.vcfanno_toml,
            params.vcfanno_resources,
            params.vcfanno_lua_func,
            params.vcf2tsv_preferred_transcript,
            params.vcf2tsv_config
        )
        ch_versions.mix(ANNOTATE.out.versions).set{ch_versions}
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )



/*


    if (publish_dir_mode == "move"){
            ch_sample_gvcf
            .map{it[0]}
            .join(ch_normal_bams)
            .mix(
                ch_sample_type.tumor
                .map{it[0]}
                .join(ch_all_bams)
            ).set{ch_output_files}
        }
        else{
            ch_all_bams.set{ch_output_files}
        }

        ch_output_files.subscribe{
            file("${params.outdir}/bams/").mkdirs()

            println("Migrating bam fils for the sample ${it[0]}")
            if (params.publish_dir_mode == "move"){
                it[1].moveTo("${params.outdir}/bams/${it[1].getName()}")
                it[2].moveTo("${params.outdir}/bams/${it[2].getName()}")
            }else{
                it[1].copyTo("${params.outdir}/bams/${it[1].getName()}")
                it[2].copyTo("${params.outdir}/bams/${it[2].getName()}")
            }
        }

*/
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
