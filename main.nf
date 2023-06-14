#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//Pipiline version
pipelineVersion = '1.0'


include {TRIMMING} from "$projectDir/modules/trimming"
include {ASSEMBLY_SHOVILL; ASSEMBLY_ASSESS; ASSEMBLY_QC } from "$projectDir/modules/assembly"
include { INDEX_REFERENCE; MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC } from "$projectDir/modules/coverage"
include { startMessage; helpMessage } from "$projectDir/modules/messages"


// Start message
startMessage(pipelineVersion)

if (params.help) {
        helpMessage()
    }else if (params.reads && params.reference && params.output ){
  workflow {

    //create raw reads channel
    raw_read_pairs_ch = Channel.fromFilePairs("$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}", checkIfExists: true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }

    //Do Trimming using Fastp
    TRIMMING(raw_read_pairs_ch)

    //Extract number of bases left after timming from fastp output i.e .json
    //READ_QC(TRIMMING.out.json, params.length_low, params.depth)

    //Do assembly and 
    ASSEMBLY_ch = ASSEMBLY_SHOVILL(TRIMMING.out.processed_reads, params.min_contig_length)
    
    
    //Assembly assess using quast
    ASSEMBLY_ASSESS(ASSEMBLY_ch)

    //Evaluate depth by the number of bases and assembly length = Assembly length/number of bases
    ASSEMBLY_QC(
        ASSEMBLY_ASSESS.out.report,
        TRIMMING.out.json,            
        params.depth
    )

    //Indexing reference using BWA
    INDEX_REFERENCE(file(params.reference))

    //Map reads using BWA MEM
    MAPPING(INDEX_REFERENCE.out, TRIMMING.out.processed_reads)

    //convert sam to bam and count sequence coverage and depth
    SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)

    //Call SNP using BCF tools and count HET SNP
    SNP_CALL(params.reference, SAM_TO_SORTED_BAM.out.bam, params.lite) | HET_SNP_COUNT
    
    //Genrate QC report 
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
    .join(ASSEMBLY_QC.out.info, failOnDuplicate: true, remainder: true)
       .map { (it[-1] == null) ? it[0..-2] + ['_'] * 3: it }
    .map { it.join',' }
    .collectFile(
        name: 'Depth_Coverage_results.csv',
        storeDir: "$params.output",
        seed: [
            'Sample_ID',
            'Mapping_QC',
            'Ref_Cov_%', 'Het-SNP#', 'Mapping_Depth', 'Assembly_Depth', 'Assembly_Depth_QC'
            ].join(','),
        sort: { it.split(',')[0] },
        newLine: true
    )
  }
} else {
    error "Please specify reads, a reference and an output directory with --reads, --reference and --output"
}
