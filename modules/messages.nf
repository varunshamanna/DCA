// Start message
void startMessage(String pipelineVersion) {
    log.info( 
        $/
        COVERGE AND DEPTH ASSESMENT PIPELINE
        |${String.format('              ║  v %-5s║', pipelineVersion)}
       /$.stripMargin()
    )
}

// Help message
void helpMessage() {
    log.info(
        '''
        |This is a Nextflow Pipeline for evaluating depth and coverage of  sequencing raw reads (FASTQ files)
        |
        |Usage:
        |
        |
        |All options are optional, some common options:
        |--reads [PATH]    Path to the input directory that contains the reads to be processed
        |--output [PATH]   Path to the output directory that save the results
        |--reference       Path to species specific reference file to be used
        |--version         To print version of this pipeline
        |
        |For all available options, please refer to README.md
        '''.stripMargin()
    )
}