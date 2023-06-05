// Start message
void startMessage(String pipelineVersion) {
    log.info( 
        $/
        ||--------------------------------------|
        || ___          ___                     |
        |||   \        /            /\          |
        |||    |      |            /__\         |
        |||__ /        \___       /    \        |
        ||                                      |
        ||DEPTH AND COVERAGE ASSESMENT PIPELINE |
        ||--------------------------------------|

        |${String.format('             ║ v %-5s║', pipelineVersion)}
       /$.stripMargin()
    )
}

// Help message
void helpMessage() {
    log.info(
        '''
        |This is a Nextflow Pipeline for evaluating depth and coverage of  sequencing raw reads (FASTQ files)
        |
        |Usage: nextflow run main.nf --reads path/to/raw/reads --reference path/to/reference/fasta --output path/to/output/directory
        |
        |
        |All options are mandatory, the options are:
        |--reads [PATH]    Path to the input directory that contains the reads to be processed
        |--output [PATH]   Path to the output directory that save the results
        |--reference       Path to species specific reference file to be used
        |--version         To print version of this pipeline
        |
        '''.stripMargin()
    )
}