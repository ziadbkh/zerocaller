
// NOTE: there is an nfcore: https://nf-co.re/modules/ENSEMBLVEP
process VEP {
    tag "$sample_id"
    label 'process_medium'

    conda "bioconda::ensembl-vep=110.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:110.0--pl5321h2a3209d_0' :
        'biocontainers/ensembl-vep:110.0--pl5321h2a3209d_0' }"

    input:
    tuple val(sample_id), path(vcf)
    //val   genome  this was removed from the script --assembly $genome \\ 
    val  species
    val   cache_version
    path  cache
    path(fasta)
    path  extra_files

    output:
    tuple val(sample_id), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(sample_id), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(sample_id), path("*.json.gz") , optional:true, emit: json
    path "*.summary.html"              , emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${sample_id}_vep"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
    
 
//&& bgzip ${prefix}.vep.vcf \
//&& tabix ${prefix}.vep.vcf.gz
