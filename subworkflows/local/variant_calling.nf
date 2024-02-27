
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { ALIGNANDSORT                  } from '../../modules/local/alignandsort'
include { MERGEBAMS                     } from '../../modules/local/mergebams.nf'
include { HAPLOTYPECALLER               } from '../../modules/local/haplotypecaller.nf'
include { MERGEVCFS                     } from '../../modules/local/mergevcfs.nf'
include { SPLITFASTQ                     } from '../../modules/local/splitfastq.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow VARIANT_CALLING {

    take:
        input
        fasta
        intervals
        fastq_max_size

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //

        Channel
            .fromPath(fasta, checkIfExists:true).first().set{ch_reference_genome}

        Channel
            .fromPath( fasta + '.{bwt,sa,ann,amb,pac}', checkIfExists:true)
            .toSortedList()
            .set{ch_reference_genome_extra_bwa}


        ch_reference_genome.map{"${it.getParent()}/${it.getBaseName()}.dict"}
        .combine(Channel
                .fromPath(fasta + '.fai', checkIfExists:true)
        )
        .set{ch_reference_genome_extra_gatk}


        Channel
            .fromPath(intervals + '*scattered.interval_list', checkIfExists:true)
            .set{ch_intervals}

        intervals_files = file(intervals + '*scattered.interval_list', checkIfExists:true)

        Channel
            .fromPath(input, checkIfExists:true )
            .splitCsv(header: true)
            .map {
                row -> {
                if (row.fastq2)
                    {
                        [row.sample_id,
                        row.type.toLowerCase(),
                        true,
                        [file(row.fastq1, checkIfExists:true), file(row.fastq2, checkIfExists:true)] ]
                    }else{
                        [row.sample_id, row.type.toLowerCase(), false, [file(row.fastq1, checkIfExists:true)] ]
                    }
                }
            }
            .set{ch_meta_complete}


        ch_meta_complete
        .map{[it[0], it[2], it[3]]}
        .branch {
            small: it[2][0].size() <= fastq_max_size
            large: it[2][0].size() > fastq_max_size
        }
        .set {
            ch_meta_all
        }

        ch_meta_complete
        .map{[it[0], it[3], Math.ceil(it[3][0].size() / fastq_max_size).intValue()]}
        .groupTuple()
        .map{[it[0], it[2].sum()]}
        .set{ch_sample_file_cnt}

        //ch_sample_file_cnt.view()
        //ch_meta_all.large.view()
        //ch_meta_all.small.view()

        ch_meta_complete
        .map{[it[0], it[1]]}
        .branch {
            tumor: it[1] == "t" || it[1] == "tumor" || it[1] == "0"
            normal: true
        }
        .set{ch_sample_type}


        SPLITFASTQ(
            ch_meta_all.large,
            Channel.value(fastq_max_size)
        )


        SPLITFASTQ.out.splited_fastq.
            transpose()
            .map{[it[0], it[2].getName().replace("R1", "").replace("R2", ""), it[1], it[2]]}
            .groupTuple(by:[0, 1, 2], size:2, remainder:true)
            .map{[it[0], it[2], it[3]]}
            .mix(ch_meta_all.small)
            .set{ch_all_fastq}

        //ch_all_fastq.view()

        ALIGNANDSORT(
            ch_all_fastq,
            ch_reference_genome.combine(ch_reference_genome_extra_bwa).flatten().toSortedList()
        )

        ch_sample_file_cnt
        .cross(ALIGNANDSORT.out.aligned_and_sorted_bam)
        .map{ tuple(groupKey(it[0][0], it[0][1]), it[1][1], it[1][2] )}
        .groupTuple()
        .branch {
            singelton: it[1].size() <= 1
            multiple: it[1].size() > 1
        }
        .set {
            ch_sample_bams
        }

        //ch_sample_bams.multiple.view()
        //ch_sample_bams.singelton.view()

        MERGEBAMS(
            ch_sample_bams.multiple
        )


        MERGEBAMS.out.merged_bam
        .mix(
            ch_sample_bams
            .singelton
            .map{[it[0], it[1][0], it[2][0]]}
        ).set{
            ch_all_bams
        }

        ch_all_bams
        .join(ch_sample_type.normal)
        .map{[it[0], it[1], it[2]]}
        .set{ch_normal_bams}


        HAPLOTYPECALLER(
            ch_normal_bams,
            ch_reference_genome
                .combine(ch_reference_genome_extra_gatk)
                .flatten()
                .toSortedList(),
            ch_intervals

        )

        HAPLOTYPECALLER.out.gvcf
        .groupTuple(size: intervals_files.size())
        .set{ch_sample_gvcf}

        MERGEVCFS(
            ch_sample_gvcf
        )


    emit:
        bams     = ch_all_bams
        bams_singelton_parts = ch_sample_bams.singelton
        bams_multiple_parts = ch_sample_bams.multiple
        gvcf     = MERGEVCFS.out.merged_vcf.join(MERGEVCFS.out.merged_vcf_tbi)
        gvcf_parts     = HAPLOTYPECALLER.out.gvcf
        versions = ch_versions
}

