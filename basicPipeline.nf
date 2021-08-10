#!/usr/bin/env nextflow
 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "/scratch/general/lustre/u0854535/RNASeqData/*.fastq.gz"
params.annot = "/scratch/general/lustre/u0854535/RNASeqRef" /*TODO: Add path to .gtf file*/
params.genome = "/scratch/general/lustre/u0854535/RNASeqRef" /*TODO: Add path to genome file*/
params.outdir = 'results'
 
/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
/*Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }
 */ 
/*
 * Step 1. Builds the genome index required by the mapping process
 */
/*
process buildIndex {
    tag "$genome.baseName"
     
    input:
    path genome from params.genome
      
    output:
    path 'genome.index*' into index_ch
        
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}
  */
/*
 * Step 2. Remove optical duplicates common in NovaSeq reads
 */
idReadFq = Channel
    .fromPath( "${params.reads}" )
    .map { file ->
        fName    = file.baseName
        id       = fName.tokenize('.')[0].tokenize('_')[0]
        read_num = fName.tokenize('.')[0].tokenize('_')[1]
        [ id, read_num, file ]

idReadFq
    .groupTuple()
    .set { reads_in }

process removeOpticalReads
    tag "$reads.baseName"

    input:
    set val(sample_id), val(reads), val(fqs) from reads_in

    output:
    path '18992X1.clump*' into reads_ch
    
    script:

    if (fqs[0].contains("_R1")) {
        fq1 = fqs[0]
        fq2 = fqs[1]
        r1  = reads[0]
        r2  = reads[1]
    } else {
        fq1 = fqs[1]
        fq2 = fqs[0]
        r1  = reads[1]
        r2  = reads[0]
    }


    """
    clumpify.sh in1=${fq1} in2=${fq2} \\
        out1=${sample_id}_${r1}.clump1.fq.gz out2=${sample_id}_${r1}.clump2.fq.gz \\
        dupedist=12000 dedupe=t optical=t 
    """

/*
 * Step 3. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
    tag "$pair_id"
      
    input:
    path genome from params.genome
    path annot from params.annot
    path index from index_ch
    tuple val(pair_id), path(reads) from read_pairs_ch
  
    output:
    set pair_id, "accepted_hits.bam" into bam_ch
  
    """
    tophat2 -p ${task.cpus} --GTF $annot genome.index $reads
    mv tophat_out/accepted_hits.bam .
    """
}
   
/*
 * Step 4. Assembles the transcript by using the "cufflinks" tool
 */
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy' 
        
    input:
    path annot from params.annot
    tuple val(pair_id), path(bam_file) from bam_ch
      
    output:
    tuple val(pair_id), path('transcript_*.gtf')
  
    """
    cufflinks --no-update-check -q -p $task.cpus -G $annot $bam_file
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}
