process ALIGNANDSORT {
    tag "${sample_id}"    
    
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'zhangb1/biobambam2-samtools-picard-bwa-samblaster-sambamba' }"
    
    
    input:
    tuple val(sample_id), val(paired_end), path(fastqs)
    path (reference_genome)
    
    output:
    tuple val(sample_id), path("${output_bam_basename}.bam"), path("${output_bam_basename}.bam.bai"), emit: aligned_and_sorted_bam
    
    
    script:
    output_bam_basename = "${sample_id}.dedup.sorted.hs38"
    metrics_filename="${sample_id}.bam.dupmetrics"
    
    fastqs_str = fastqs.join(" ")
    bamsorThreads = Math.min(12,task.cpus)
    refernce_genome_main = reference_genome.find { it.getName().endsWith(".fa") }

    """
    bwa mem \
        -t ${task.cpus} \
        -R @RG\\\\tID:"${sample_id}"\\\\tLB:None\\\\tPL:ILLUMINA\\\\tPU:NA\\\\tSM:"${sample_id}" \
        ${task.ext.bwa_args} \
        ${refernce_genome_main} \
        ${fastqs_str} | \
    bamsormadup \
        threads=${bamsorThreads} \
        ${task.ext.bamsormadup_args} \
        reference=${refernce_genome_main} \
        tmpfile=tmp_${output_bam_basename} \
        fragmergepar=${bamsorThreads} \
        M=${metrics_filename} \
        indexfilename=${output_bam_basename}.bam.bai \
        > ${output_bam_basename}.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        bamsormadup: \$(echo \$(bamsormadup --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
    END_VERSIONS
    """

}