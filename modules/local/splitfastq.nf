
process SPLITFASTQ {
    
    conda "bioconda::fastqsplitter=1.2.0 "
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqsplitter:1.2.0--py39hbf8eff0_3' :
        'biocontainers/fastqsplitter:1.2.0--py310h1425a21_3' }"
        
    input:
        tuple val(sample_id), val(paired_end), path(fastqs)
        val(chunk_size)
    output: 
        tuple val(sample_id), val(paired_end), path("*.fastq.gz"), emit: splited_fastq
    

    script:
    fastq1 = fastqs[0]
    chuncks_cnt = Math.ceil(fastq1.size() / chunk_size).intValue()

    read1_out = (1..chuncks_cnt).collect{ "-o " + fastq1.getSimpleName() + "." + it + ".fastq.gz"}
        
    
    if (paired_end){
        fastq2 = fastqs[1]
        read2_out = (1..chuncks_cnt).collect{"-o " + fastq2.getSimpleName() + "." + it + ".fastq.gz"}
        read2_command = "fastqsplitter -i ${fastq2} -t ${task.cpus} ${task.ext.args} ${read2_out.join(' ')}"
    }else{
        read2_command = ""
    }

    """
    fastqsplitter -i ${fastq1} -t ${task.cpus} ${task.ext.args}  ${read1_out.join(' ')}
    ${read2_command}
    
    """
}