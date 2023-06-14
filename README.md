# Depth & Coverage Assesment Pipeline <!-- omit in toc -->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-23.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/singularity/)


This pieline can be used to evaluate Depth, Coverage and Heterozygous SNP count to evaluate the mapping QC of raw reads.

Requirement:
Nextflow
Docker

## Accepted Inputs
- Currently, only Illumina paired-end short reads are supported
- Each sample is expected to be a pair of raw reads following this file name pattern: 
  - `*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}` 
    - example 1: `SampleName_R1_001.fastq.gz`, `SampleName_R2_001.fastq.gz`
    - example 2: `SampleName_1.fastq.gz`, `SampleName_2.fastq.gz`
    - example 3: `SampleName_R1.fq`, `SampleName_R2.fq`

## Setup 
1. Clone the repository (if Git is installed on your system)
    ```
    git clone https://github.com/varunshamanna/DCA.git
    ```
    or 
    
    Download and unzip the [repository](https://github.com/varunshamanna/DCA/archive/refs/heads/master.zip)
2. Go into the local copy of the repository
    ```
    cd DCA
    ```
## Run
> ⚠️ Docker  must be running.
<!-- -->
> ⚠️ If this is the first run , an Internet connection is required.
<!-- -->
> ℹ️ By default, Docker is used as the container engine and all the processes are executed by the local machine. See [Profile](#profile) 

The pipeline requires three inputs which are mandatory

Path to the Illumina paried-end Raw reads (The folder must have raw reads of specific species)
  ```
  --reads /path/to/raw-reads-directory
  ```

A standard species specific reference fasta file (specific to your raw reads)
  ```
  --reference /path/to/reference/fasta
  ```

Path to a output folder to save the outputs
  ```
  --output path/to/output/
  ```

Running the pipeline

  ```
  nextflow run main.nf --reads /path/to/raw/reads-directory --reference /path/to/reference/fasta --output path/to/output -resume
  ```
 
Example: 

  ```
  nextflow run main --reads /data/staph_aureus  --reference /data/references/SAU/SAU_reference.fna --output sau_dca -resume
  ```

# Credits
The concept for developing this pipeline was deveoped from [GPS-Unified-Pipeline](https://github.com/HarryHung/gps-unified-pipeline) built by [@harry_hung](https://github.com/HarryHung) for the [GPS project](https://www.pneumogen.net/gps/).


