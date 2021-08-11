#!/usr/bin/env nextflow
 
params.reads = "/uufs/chpc.utah.edu/common/home/u0854535/practice/nextflow-tutorial/data/ggal/gut_{1,2}.fq"
params.outDir = "/uufs/chpc.utah.edu/common/home/u0854535/workflows/RNASeq/RNASeq_Pipeline/out"
params.genomeRef = "/uufs/chpc.utah.edu/common/home/u0854535/RNASeq/scripts/star125"

Channel.fromFilePairs( "${params.reads}", checkIfExists:true )
    .into { reads_in; reads_in2; reads_in3 }
    
reads_in.view()
/**********Remove optical reads if using HiSeq 2500/MiSeq/NextSeq data**********
process removeOpticalReads {

    publishDir params.outDir

    input:
        tuple val(sample_id), path(fqs) from reads_in2

    output:
        tuple path("${sample_id}_1.clump1.fq.gz"), path("${sample_id}_2.clump2.fq.gz") into reads_ch
    
    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        println "reads1 $fq1\nreads2 $fq2\nsampleID $sample_id"

        """
        clumpify.sh in1=${fq1} in2=${fq2} out1=${sample_id}_1.clump1.fq.gz out2=${sample_id}_2.clump2.fq.gz dupedist=12000 dedupe=t optical=t 
        """
}
*/

process fastqc {

    publishDir params.outDir    

    input:
        tuple val(sample_id), path(fqs) from reads_in2

    output:
        file '*.html' into qc_ch

    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        println "reads1 $fq1\nreads2 $fq2\nsampleID $sample_id"
        
        """
        module load fastqc
        fastqc -f fastq $fq1 > result1
        fastqc -f fastq $fq2 > result2
        """
}

process sortFastqReads {

    publishDir params.outDir

    input:
        tuple val(sample_id), path(fqs) from reads_in3

    output:
        /*file '*.fastq' into sorted_reads_ch*/
        tuple val(sample_id), path('*.fastq') into sorted_reads_ch

    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        """
        cat $fq1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sample_id}sorted_1.fastq
        cat $fq2 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sample_id}sorted_2.fastq
        """
}

process alignReads {

    publishDir params.outDir

    input:
        /*tuple val(sample_id), path(fqs) from reads_in3*/
        tuple val(sample_id), path(fqs) from sorted_reads_ch

    output:
        file '*.bam' into aligned_ch 
    
    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        """
        module load gcc/8.3.0
        STAR --genomeDir $params.genomeRef \
             --runThreadN 24 \
             --readFilesIn $fq1 $fq2 \
             --twopassMode Basic \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM \
             --outWigType bedGraph \
             --outWigStrand Unstranded
        """
}
