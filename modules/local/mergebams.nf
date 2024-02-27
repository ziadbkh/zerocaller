
process MERGEBAMS {
    tag "${sample_id}"

    conda  "bioconda::sambamba=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0--h98b6b92_0':
        'biocontainers/sambamba:1.0--h98b6b92_0' }"
        
    input:
    tuple val(sample_id), path (bam_list, stageAs: 'file??.bam'), path(bam_bai_list, stageAs: 'file??.bam.bai')
    
    output:
    tuple val(sample_id), path ("${output_bam}"), path ("${output_bam}.bai"), emit: merged_bam
    
    script:
    output_bam="${sample_id}.dedup.merged.sorted.hs38.bam"
    
    """
    sambamba merge -t ${task.cpus} ${output_bam} ${bam_list.join(" ")} ${task.ext.args}
    """
}