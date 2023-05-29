// Run fastp to preprocess the FASTQs
process TRIMMING {
    label 'fastp_container'

    tag "$sample_id"

    publishDir "${params.output}/fastp_out", mode: 'link', pattern: "*.html"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(processed_one), path(processed_two), path(processed_unpaired), emit: processed_reads
    tuple val(sample_id), path(html), emit: html

    script:
    read_one="${reads[0]}"
    read_two="${reads[1]}"
    processed_one="processed-${sample_id}_1.fastq.gz"
    processed_two="processed-${sample_id}_2.fastq.gz"
    processed_unpaired="processed-${sample_id}_unpaired.fastq.gz"
    html="${sample_id}.fastp.html"
    """
    fastp --thread 8 --in1 "$read_one" --in2 "$read_two" --out1 "$processed_one" --out2 "$processed_two" --unpaired1 "$processed_unpaired" --unpaired2 "$processed_unpaired"
    mv *.html $html
    """
}