
process VCFANNO {
    tag "$sample_id"
    label 'process_low'

    conda "bioconda::vcfanno=0.3.3"
    /*
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfanno:0.3.3--h9ee0642_0':
        'biocontainers/vcfanno:0.3.3--h9ee0642_0' }"
    */
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcfanno:0.3.3--h9ee0642_0':
        'pgc-images.sbgenomics.com/syan/vcfanno:production' }"
    
    input:
    tuple val(sample_id), path(vcf)
    path toml
    path lua
    path resources

    output:
    tuple val(sample_id), path("*.vcf.gz") , emit: vcf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}-vcfanno"
    def lua_cmd = lua ? "--lua ${lua}" : ""
    """
    ${args2} vcfanno \\
        -p ${task.cpus} \\
        ${args} \\
        ${lua_cmd} \\
        ${toml} \\
        ${vcf} \\
        > ${prefix}.vcf

    bgzip --threads ${task.cpus} ${prefix}.vcf
    tabix ${prefix}.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfanno: \$(echo \$(vcfanno 2>&1 | grep version | cut -f3 -d' ' ))
    END_VERSIONS
    """
}