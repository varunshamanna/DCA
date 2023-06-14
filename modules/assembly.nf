// Return sample_id and assembly, and hardlink the assembly to ${params.output}/assemblies directory
process ASSEMBLY_SHOVILL {
    label 'shovill_container'
 
    tag "$sample_id"

    publishDir "${params.output}/assemblies", mode: 'link'

    input:
    tuple val(sample_id), path(read1), path(read2), path(unpaired)
    val min_contig_length

    output:
    tuple val(sample_id), path(fasta)

    script:
    fasta="${sample_id}.contigs.fasta"
    """
    shovill --R1 "$read1" --R2 "$read2" --outdir results --cpus 8 --minlen "$min_contig_length" --force
    mv results/contigs.fa "${fasta}"
    """
}

// Run quast to assess assembly quality
process ASSEMBLY_ASSESS {
    label 'quast_container'

    tag "$sample_id"

    publishDir "${params.output}/quast_reports", mode: 'link'

    input:
    tuple val(sample_id), path(assembly)
    output:
    tuple val(sample_id), path('results/*_quast_report.tsv'), emit: report

    script:
    """
    quast.py -o results "$assembly"
    mv results/report.tsv results/"$sample_id"_quast_report.tsv
    """
}

// Extract assembly QC information and determine QC result based on report.tsv from Quast, and total base count
process ASSEMBLY_QC {
    label 'bash_container'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(report)
    tuple val(sample_id), path(json)
    val(qc_depth)

    output:
    tuple val(sample_id), env (ASSEMBLY_DEPTH), env(ASSEMBLY_DEPTH_QC), emit: info

    script:
    """
    REPORT="$report"
    JSON="$json"
    QC_DEPTH="$qc_depth"
    
    source assembly_qc.sh      
    """
}