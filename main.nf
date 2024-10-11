
nextflow.enable.dsl = 2

process COPY_FASTQ{
    container "ubuntu:24.10" 

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fastq"), emit: fastq_copied
    path("abc.txt")
    path("head.txt")

    script:
    """
    touch abc.txt
    cp ${reads} ${sample_id}.fastq
    touch head.txt
    head ${reads} > head.txt
    ls -Lla
    """
}


process COPY_FASTA{
    container "ubuntu:24.10"
    
    input:
    tuple val(sample_id), path(reads)
    path(genome)

    output:
    tuple val(sample_id), path("${sample_id}*1234.f*q")
    path "${sample_id}.fasta"
    
    script:
    """
    mv ${reads} ${sample_id}1234.fq
    touch head.txt
    head ${sample_id}1234.fq > head.txt
    cp ${genome} ${sample_id}.fasta
    """
}

workflow {
    
        read_pairs_ch = Channel
            .fromPath( params.csv_input )
            .splitCsv(header: true, sep: ',')
            .map {row -> tuple(row.sample, [row.path_r1])}
            .view()

    COPY_FASTQ( read_pairs_ch )
    COPY_FASTA( COPY_FASTQ.out.fastq_copied, params.fasta)
}
