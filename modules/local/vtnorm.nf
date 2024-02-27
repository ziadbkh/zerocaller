
process VTNORM{
    tag "$sample_id"
    
    conda "bioconda::vt=2015.11.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:2015.11.10--h5ef6573_4':
        'pgc-images.sbgenomics.com/syan/vt-norm:production'}"
        //'biocontainers/vt:2015.11.10--h5ef6573_4' }"

    input:
    tuple val(sample_id), path(vcf), path(tbi) //, path(intervals)
    path(fasta)
    path(fai)

    output:
    tuple val(sample_id), path("*.vcf.gz")       , emit: vcf
    tuple val(sample_id), path("${fasta}.fai")   , emit: fai, optional: true
    path "versions.yml"                     , emit: versions

    //when: task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    //def regions = intervals ? "-i ${intervals}" : ""
    def regions = ""
    if ("$vcf" == "${prefix}.vcf" || "$vcf" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def VERSION = "2015.11.10" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vt normalize \\
        -o ${prefix}.vcf \\
        -r ${fasta} \\
        ${regions} \\
        ${args} \\
        ${vcf}

    bgzip ${args2} --threads ${task.cpus} ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt: ${VERSION}
    END_VERSIONS
    """
// lsblk
//       mount
//       df -h
//       #vt normalize -o norm.vcf -n -r /sbgenomics/Projects/647fff0e-2f15-462a-afab-0f41df0069ba/Assets/genome/hg38/GRCh38_no_alt.fa /sbgenomics/Projects/647fff0e-2f15-462a-afab-0f41df0069ba/P006105_ASXL2.vcf.gz
//       vt normalize -o output.vcf -n -r ${refgnome} ${inputvcf}
//       #mv norm.vcf P006105_ASXL2.vcf
//       #bgzip -f P006105_ASXL2.vcf
//       #tabix -f P006105_ASXL2.vcf.gz
}