// Index reference
process INDEX_REFERENCE {
  label 'bwa_container'

  input:
    file(reference)
  output:
    path("index")

  script:
    """
    mkdir index
    bwa index ${reference} -p index/reference
    """

}


// Map the reads to reference using BWA-MEM algorithm
// Return mapped SAM
process MAPPING {
    label 'bwa_container'
    tag "$sample_id"
     
    input:
    path (index)
    tuple (val(sample_id), path(read1), path(read2), path(unpaired))

    output:
    tuple val(sample_id), path(sam), emit: sam

    script:
    sam="${sample_id}_mapped.sam"
    """
    bwa mem -t 8 -a "${index}/reference" <(zcat -f -- < "$read1") <(zcat -f -- < "$read2") > "$sam"
    """
}

// Convert mapped SAM into BAM and sort it
// Return mapped and sorted BAM, and reference coverage percentage by the reads
process SAM_TO_SORTED_BAM {
    label 'samtools_container'

    tag "$sample_id"
    publishDir "${params.output}/sorted_bam", mode: 'link'

    input:
    tuple val(sample_id), path(sam)
    val lite

    output:
    tuple val(sample_id), path(bam), emit: bam
    tuple val(sample_id), env(COVERAGE), emit: ref_coverage
    tuple val(sample_id), env(DEPTH), emit: ref_depth

    script:
    bam="${sample_id}_mapped_sorted.bam"
    """
    samtools view -@ 8 -b "$sam" > mapped.bam

    samtools sort -@ 8 -o "$bam" mapped.bam
    rm mapped.bam

    if [ $lite = true ]; then
        rm `readlink -f "$sam"`
    fi

    BAM="$bam"
    source get_ref_coverage.sh
    """
}

// Return .vcf by calling the SNPs
process SNP_CALL {
    label 'bcftools_container'

    tag "$sample_id"
    publishDir "${params.output}/vcf", mode: 'link'

    input:
    path reference
    tuple val(sample_id), path(bam)
    val lite

    output:
    tuple val(sample_id), path(vcf), emit: vcf

    script:
    vcf="${sample_id}.vcf"
    """
    bcftools mpileup --threads 8 -f "$reference" "$bam" | bcftools call --threads 8 -mv -O v -o "$vcf"

    if [ $lite = true ]; then
        rm `readlink -f "$bam"`
    fi
    """
}

// Return non-cluster heterozygous SNP (Het-SNP) site count
process HET_SNP_COUNT {
    label 'python_container'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), env(OUTPUT), emit: result

    script:
    """
    OUTPUT=`het_snp_count.py "$vcf" 50`
    """
}

// Extract mapping QC information and determine QC result based on reference coverage and count of Het-SNP sites
process MAPPING_QC {
    label 'bash_container'

    tag "$sample_id"

    input:
    tuple val(sample_id), val(ref_coverage), val(ref_depth), val(het_snp_count)
    val(qc_ref_coverage)
    val(qc_het_snp_site)
    val(qc_depth)

    output:
    tuple val(sample_id), env(COVERAGE), env(HET_SNP), env(DEPTH), emit: info
    tuple val(sample_id), env(MAPPING_QC), emit: result

    script:
    """
    COVERAGE="$ref_coverage"
    HET_SNP="$het_snp_count"
    DEPTH="$ref_depth"
    QC_REF_COVERAGE="$qc_ref_coverage"
    QC_HET_SNP_SITE="$qc_het_snp_site"
    QC_DEPTH="$qc_depth"

    source mapping_qc.sh
    """
}