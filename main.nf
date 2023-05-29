#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//Pipiline version
pipelineVersion = '1.0'

//params.coverage_threshold = 5

// Start message
//startMessage(pipelineVersion)

include {TRIMMING} from "$projectDir/modules/trimming"
include { INDEX_REFERENCE; MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC } from "$projectDir/modules/coverage"
include { startMessage; helpMessage } from "$projectDir/modules/messages"


// Start message
startMessage(pipelineVersion)

if (params.help) {
        helpMessage()
    }else if (params.reads && params.reference && params.output ){
  workflow {
    raw_read_pairs_ch = Channel.fromFilePairs("$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}", checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    TRIMMING(raw_read_pairs_ch)
    INDEX_REFERENCE(file(params.reference))
    MAPPING(INDEX_REFERENCE.out, TRIMMING.out.processed_reads)
    SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)
    SNP_CALL(params.reference, SAM_TO_SORTED_BAM.out.bam, params.lite) | HET_SNP_COUNT
    MAPPING_QC(
        SAM_TO_SORTED_BAM.out.ref_coverage
        .join(SAM_TO_SORTED_BAM.out.ref_depth, failOnDuplicate: true, failOnMismatch: true)
        .join(HET_SNP_COUNT.out.result, failOnDuplicate: true, failOnMismatch: true), 
        params.ref_coverage,
        params.het_snp_site,
        params.depth)
    MAPPING_QC.out.result
    .join(MAPPING_QC.out.info, failOnDuplicate: true, remainder: true)
        .map { (it[-1] == null) ? it[0..-2] + ['_'] * 3 : it }
    .map { it.join',' }
    .collectFile(
        name: 'Depth_Coverage_results.csv',
        storeDir: "$params.output",
        seed: [
            'Sample_ID',
            'Mapping_QC',
            'Ref_Cov_%', 'Het-SNP#', 'Seq_Depth' 
            ].join(','),
        sort: { it.split(',')[0] },
        newLine: true
    )
    //ASSESS_COVERAGE(reads, INDEX_REFERENCE.out, params.coverage_threshold)
    //COMBINE_COVERAGE_RESULTS(ASSESS_COVERAGE.out.collect())
  }
} else {
    error "Please specify reads, a reference and an output directory with --reads, --reference and --output_dir"
}
