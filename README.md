# RNASeqPipeline
This is a simple pipeline for RNASeq data analysis based off of the workflow used by the Bioinformatics Core at the University of Utah HCRI. Currently input files can only be in fastq format.

While changes are frequently being made and the pipeline is being tested, to avoid using up our CHPC allocation I am running the pipeline out of the kingspeak directory. There is a `nextflow.config` file in the `/kingspeak/configs` directory that can be updated for each unique run of the pipeline however, I am keeping things simple while testing and have the parameters listed in the top of the `main.nf` file.

## Usage

`nextflow run main.nf`
