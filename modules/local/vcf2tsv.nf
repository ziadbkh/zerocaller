//Process
process VCF2TSV {    
    tag "$sample_id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'ziadbkh/vcf2tsv-sp:v1.4.1' }"

    input:
        tuple val(sample_id), path(vcf)
        path pythonconfig
        path preftsv


    output:
        tuple val(sample_id), path("*.tsv") , emit: tsv
        path "versions.yml"                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vcf2tsv --nprocs ${task.cpus} --config $pythonconfig  --prefer $preftsv  ${args} $vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2tsv: 1.4.1
    END_VERSIONS
    """
}
