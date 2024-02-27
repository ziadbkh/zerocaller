process TABIX {
    tag "$sample_id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(sample_id), path(tab)

    output:
    tuple val(sample_id), path(tab), path("*.tbi"), optional:true, emit: tbi
    tuple val(sample_id), path(tab), path("*.csi"), optional:true, emit: csi
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    tabix $args $tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}