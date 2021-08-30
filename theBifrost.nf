#!/usr/bin/env nextflow

//Check for input reads 
if (params.reads) { Channel.fromFilePairs( "${params.reads}", checkIfExists:true )
    .into { reads_in; reads_in2; reads_in3 }
} else { exit 1, 'Input reads not specified!' }

// TODO: some kind of if else that depends on whether or not a path to single end reads
// is given. It'd be cool if that could be added here to the config or over CLI as --singleEndReads

//Channel.fromFilePairs( "${params.reads}", checkIfExists:true )
//    .into { reads_in; reads_in2; reads_in3 }
    
//reads_in.view()


//TODO: Test removal of optical reads with our RNASeq data from the core
/**********Remove optical reads if using HiSeq 2500/MiSeq/NextSeq data**********
process removeOpticalReads {

    publishDir params.outDir

    input:
        tuple val(sample_id), path(fqs) from reads_in

    output:
        tuple val(sample_id), path("*.fq.gz") into reads_ch
    
    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        println "reads1 $fq1\nreads2 $fq2\nsampleID $sample_id"

        """
        clumpify.sh in1=${fq1} in2=${fq2} out1=${sample_id}_1.clump1.fq.gz out2=${sample_id}_2.clump2.fq.gz dupedist=12000 dedupe=t optical=t 
        """
}


//Trim adapters
process trimAdapt {

    input:
        tuple val(sample_id), path(fqs) from reads_ch

    output:
        tuple val(sample_id), path(*.fq.gz) into fastqs_ch

    script:
        fq1 = fqs[0]
        fq2 = fqs[1]

        """
        cutadapt -j 24 -O 6 -m 20 \
            -a CTGTCTCTTATACACATCT \
            -A CTGTCTCTTATACACATCT \
            -o ${sample_id}.1.fq.gz -p ${sample_id}.2.fq.gz \
            $fq1 $fq2
        """
}

*******************************************************************************/


//Run reads through fastqc, find out if adapters are present in results
process fast_qc {

    tag "$sample_id"

    //publishDir params.outDir    

    input:
        tuple val(sample_id), path(fqs) from reads_in2
        //The tuple looks like [sample_id, [fastq1, fastq2]]

    //output:
        //file '*.html' into qc_ch

    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        println "reads1 $fq1\nreads2 $fq2\nsampleID $sample_id"
        
        """
        module load fastqc
        fastqc -f fastq $fq1 -o $params.fastqcOut
        fastqc -f fastq $fq2 -o $params.fastqcOut
        """
}


//Sort fastq reads by sequence identifier
/*process sortFastqReads {

    tag "$sample_id"

    publishDir params.outDir

    input: //TODO: have this read in from trimmed reads out channel once working
        tuple val(sample_id), path(fqs) from reads_in3

    output:
        tuple val(sample_id), path('*.fastq.gz') into sorted_reads_ch

    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        """
        cat $fq1 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sample_id}sorted_1.fastq
        gzip ${sample_id}sorted_1.fastq
        cat $fq2 | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ${sample_id}sorted_2.fastq
        gzip ${sample_id}sorted_2.fastq
        """
}*/


//Align reads
process alignReads {

    tag "$sample_id"

    publishDir params.bams

    input:
        //tuple val(sample_id), path(fqs) from sorted_reads_ch
        tuple val(sample_id), path(fqs) from reads_in3

    //output:
        //tuple val(sample_id), path("*.bam") into aligned_ch 
    
    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        //STAR cannot read files that are compressed unless you provide the --readFilesCommand zcat
        """
        module load gcc/8.3.0
        STAR --genomeDir $params.genomeRef \
             --runThreadN 24 \
             --readFilesIn $fq1 $fq2 \
             --twopassMode Basic \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM \
             --outWigType bedGraph \
             --outWigStrand Unstranded \
             --outFileNamePrefix ${params.bams}/${sample_id} \
             --outTmpKeep All \
             --outStd Log \
             --readFilesCommand zcat
        """
}

//Set the bams to a new channel for checking alignment
Channel.from("${params.bams}/*.bam").set{ aligned_ch2 }
//TODO: Need indexes of Bams?

//Count uniquely aligned reads overlapping features in the GTF file
/*process countFeatures {

    publishDir params.outDir

    beforeScript "echo ${alignedReads}"

    input:
        tuple val(sample_id), path(alignedReads) from aligned_ch

    //output:
        //file '*.counts' into counts_ch

    script:
        """
        featureCounts -T 24 -p -s 2  --largestOverlap -a ${params.refGTF} -o ${sample_id}.counts $alignedReads
        """
}*/


//Check alignment QC
/*process checkAlignment {

    input:
        file(inputBam) from aligned_ch2

    script:
        """
        java -jar /uufs/chpc.utah.edu/common/home/u0854535/software/picard-2.25.6/picard/build/libs CollectRnaSeqMetrics \
            I=$inputBam \
            O=output.RNA_Metrics \
            REF_FLAT=$params.refFlat \
            STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
        """ 
        //Do we want to include ribsomal intervals in this check?
        

}


process indexStats {

    script:
        """
        samtools idxstats
        """
}*/
//R/4.0.2 has DESeq2 installed in its site library (.Library.site) 
