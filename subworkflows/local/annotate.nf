/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

//include {VCF2TSV} from '../modules/local/vcf2tsv.nf'

include {VTNORM} from '../../modules/local/vtnorm.nf'
include {VEP} from '../../modules/local/vep.nf'
include {VCFANNO} from '../../modules/local/vcfanno.nf'
include {VCF2TSV} from '../../modules/local/vcf2tsv.nf'
include {TABIX} from '../../modules/local/tabix.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow ANNOTATE {

    take:
        vcfs
        fasta
        vep_species
        vep_cache_version
        vep_cache
        vep_plugins
        vcfanno_toml
        vcfanno_resources
        vcfanno_lua_func
        vcf2tsv_preferred_transcript
        vcf2tsv_config
    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //

        Channel
            .fromPath( fasta, checkIfExists:true).first().set{ch_reference_genome}

        Channel
            .fromPath(fasta + '.{bwt,sa,ann,amb,pac}', checkIfExists:true).toSortedList().set{ch_reference_genome_extra_bwa}


        ch_reference_genome.map{file("${it.getParent()}/${it.getBaseName()}.dict", checkIfExists:true)}
        .combine(Channel
                .fromPath( fasta + '.fai', checkIfExists:true)
        ).first()
        .set{ch_reference_genome_extra_gatk}

        Channel.fromPath(vep_cache, checkIfExists:true).set{ch_vep_cache}
        Channel.fromPath(vep_plugins, checkIfExists:true).set{ch_vep_plugins}
        Channel.fromPath(vcfanno_toml, checkIfExists:true).first().set{ch_vcfanno_toml}
        Channel.fromPath("${vcfanno_resources}/*", checkIfExists:true).toSortedList().set{ch_vcfanno_resources}
        Channel.fromPath(vcfanno_lua_func, checkIfExists:true).first().set{ch_vcfanno_lua_func}
        Channel.fromPath(vcf2tsv_preferred_transcript, checkIfExists:true).first().set{ch_vcf2tsv_preferred_transcript}
        Channel.fromPath(vcf2tsv_config, checkIfExists:true).first().set{ch_vcf2tsv_config}

        TABIX(
            vcfs.filter{it.size() == 2}
        )

        vcfs.view()
        //ch_reference_genome.view()
        //Channel.fromPath( fasta + '.fai', checkIfExists:true).view()
        VTNORM(
            vcfs.filter{it.size() == 3}
            .mix(TABIX.out.tbi),
            ch_reference_genome,
            Channel.fromPath( fasta + '.amb', checkIfExists:true).first()
        )

        VEP(
            VTNORM.out.vcf,
            Channel.value(vep_species),
            Channel.value(vep_cache_version),
            ch_vep_cache.first(),
            ch_reference_genome,
            ch_vep_plugins.first()
        )

        VCFANNO(
            VEP.out.vcf,
            ch_vcfanno_toml,
            ch_vcfanno_lua_func,
            ch_vcfanno_resources
        )

        VCF2TSV(
            VCFANNO.out.vcf,
            ch_vcf2tsv_config,
            ch_vcf2tsv_preferred_transcript
        )

    emit:
        tsvs= VCF2TSV.out.tsv
        vcf_normalised= VTNORM.out.vcf
        vcf_vep = VEP.out.vcf
        vcf_annotated = VCFANNO.out.vcf
        versions = ch_versions

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
