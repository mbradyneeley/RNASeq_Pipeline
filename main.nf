#!/usr/bin/env nextflow

   // params.pairedEndReads   = "/scratch/general/lustre/u0854535/RNASeqData/*_{1,2}.fastq.gz"
   // //params.pairedEndReads   = "/scratch/general/lustre/u0854535/RNASeqData/JC7904_{1,2}.fastq.gz"
   // params.outDir       = "/uufs/chpc.utah.edu/common/home/u0854535/workflows/RNASeq/RNASeq_Pipeline/kingspeak/out"
   // params.genomeRef    = "/uufs/chpc.utah.edu/common/home/u0854535/RNASeq/scripts/star125"
   // params.refGTF       = "/scratch/general/lustre/u0854535/RNASeqRef/Homo_sapiens.GRCh38.104.gtf"
   // params.refFA        = "/uufs/chpc.utah.edu/common/home/u0854535/referenceFiles/hg38_v0_Homo_sapiens_assembly38.fasta"
   // 
   // //Flat file for hg38/GR38
   // params.refFlat      = "/scratch/general/lustre/u0854535/RNASeqRef/refFlat_hg38.txt"
   // 
   // //Make sure the scratch directory is up to date
   // params.fastqcOut    = "/scratch/general/lustre/u0854535/theBifrost/results/fastqc"
   // params.bams         = "/scratch/general/lustre/u0854535/theBifrost/results/bam"
   // params.multi_qc     = "/scratch/general/lustre/u0854535/theBifrost/results/multi_qc"
   // 
   // //DESeq2 parameters
   // params.DESeq        = "/scratch/general/lustre/u0854535/theBifrost/results/deseq2"
   // params.phenoFile    = "/uufs/chpc.utah.edu/common/home/u0854535/workflows/RNASeq/RNASeq_Pipeline/kingspeak/RNAseq_submission_72samples_pheno.xlsx"
   // 
   // params.cpusLeft     = 5
   // 
   // //Project directory for accessing local files
   // params.projectDir   = "/uufs/chpc.utah.edu/common/home/u0854535/workflows/RNASeq/RNASeq_Pipeline/kingspeak"
    
log.info """\
===========================================================================

██████  ███    ██  █████  ███████ ███████  ██████  
██   ██ ████   ██ ██   ██ ██      ██      ██    ██ 
██████  ██ ██  ██ ███████ ███████ █████   ██    ██ 
██   ██ ██  ██ ██ ██   ██      ██ ██      ██ ▄▄ ██ 
██   ██ ██   ████ ██   ██ ███████ ███████  ██████  
                                              ▀▀                                                               

██████  ██ ██████  ███████ ██      ██ ███    ██ ███████ 
██   ██ ██ ██   ██ ██      ██      ██ ████   ██ ██      
██████  ██ ██████  █████   ██      ██ ██ ██  ██ █████   
██      ██ ██      ██      ██      ██ ██  ██ ██ ██      
██      ██ ██      ███████ ███████ ██ ██   ████ ███████ 

P E Z Z O L E S I   L A B   2 0 2 1
B R A D Y   N E E L E Y
===========================================================================

    reads        = $params.pairedEndReads
    outDir       = $params.outDir
    genomeRef    = $params.genomeRef
    refGTF       = $params.refGTF
    refFlat      = $params.refFlat
    fastqcOut    = $params.fastqcOut
    bams         = $params.bams
    multi_qc     = $params.multi_qc
    DESeq        = $params.DESeq
    sampleIDs    = $params.phenoFile

==========================================================================
"""                

//Check for input reads 
if (params.pairedEndReads) { Channel.fromFilePairs( "${params.pairedEndReads}", checkIfExists:true )
    .into { reads_in; reads_in2; reads_in3 }
} else { exit 1, 'Input reads not specified!' }

//Channel.fromFilePairs( "${params.pairedEndReads}", checkIfExists:true )
//    .into { reads_in; reads_in2; reads_in3 }
    
//reads_in.view()


//TODO: Test removal of optical reads with our RNASeq data from the core
process removeOpticalReads {

    tag "$sample_id"

    publishDir params.outDir

    input:
        tuple val(sample_id), path(fqs) from reads_in

    output:
        tuple val(sample_id), path("*.fq.gz") into reads_ch
    
    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        //println "reads1 $fq1\nreads2 $fq2\nsampleID $sample_id"

        """
        clumpify.sh in1=${fq1} in2=${fq2} out1=${sample_id}_1.clump1.fq.gz out2=${sample_id}_2.clump2.fq.gz dupedist=12000 dedupe=t optical=t 
        """
}


//Trim adapters
process trimAdapt {

    tag "$sample_id"

    input:
        tuple val(sample_id), path(fqs) from reads_ch

    output:
        tuple val(sample_id), path("*.fq.gz") into fastqs_ch

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

fastqs_ch.into { trimmedFQs1; trimmedFQs2 }


//Run reads through fastqc, find out if adapters are present in results
process fast_qc {

    scratch '/scratch/general/lustre/u0854535/temp'

    tag "$sample_id"

    //publishDir params.outDir    

    input:
        tuple val(sample_id), path(fqs) from trimmedFQs1
        //tuple val(sample_id), path(fqs) from reads_in2
        //The tuple looks like [sample_id, [fastq1, fastq2]]

    output:
        val true into done_fastqc
        //file '*.html' into qc_ch

    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        //println "reads1 $fq1\nreads2 $fq2\nsampleID $sample_id"
        
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

    scratch '/scratch/general/lustre/u0854535/temp'

    tag "$sample_id"

    publishDir "${params.bams}", mode: 'copy', pattern: '*.bam'
    //publishDir params.bams

    input:
        //tuple val(sample_id), path(fqs) from sorted_reads_ch
        //tuple val(sample_id), path(fqs) from reads_in3
        tuple val(sample_id), path(fqs) from trimmedFQs2
        val flag from done_fastqc

    output:
        path("*.sortedByCoord.out.bam") into aligned_ch 
    
    script:
        fq1 = fqs[0]
        fq2 = fqs[1]
        //STAR cannot read files that are compressed unless you provide the --readFilesCommand zcat
        //limitBAMsortRAM was 50000000000 or 50GB now its 140GB
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
             --outFileNamePrefix ${sample_id} \
             --outTmpKeep All \
             --outStd Log \
             --limitBAMsortRAM 140000000000 \
             --readFilesCommand zcat
        """
}


process sortBams {

    publishDir "${params.bams}", mode: 'copy', pattern: '*.bam'

    input:
        path(bam) from aligned_ch

    output:
        path(bam) into sortedAligned_ch 
    
    script:
        """
        ml samtools
        samtools sort $bam -o ${bam.simpleName}_sorted.bam
        """

}

//Set the bams to a new channel for checking alignment
//Channel.fromPath("${params.bams}/*sortedByCoord.out.bam").into{ aligned_ch3; aligned_ch2; aligned_ch }
sortedAligned_ch.into{ aligned_ch4;aligned_ch3; aligned_ch2; aligned_ch1 }
//TODO: Need indexes of Bams?
//aligned_ch2.view()
//Count uniquely aligned reads overlapping features in the GTF file
//JUST PASS ALL THE BAMS INTO FEATURECOUNTS AT ONCE, THEN FEEDING TO DESEQ2 IS EASY
process countFeatures {

    publishDir params.outDir

    beforeScript "echo ${alignedReads}"

    input:
        path '*.bam' from aligned_ch1.toList()

    output:
        //this works file '*.counts.txt' into counts_ch
        path '*.counts.txt' into counts_ch1

    script:
        //println "$bam"
        """ 
        featureCounts -T 24 -p -s 2  --largestOverlap -a ${params.refGTF} -o featCounts.counts.txt *.bam
        """
        //featureCounts -T 24 -p -s 2  --largestOverlap -a ${params.refGTF} -o ${bam.simpleName}.counts.txt $bam
}

counts_ch1.into{ counts_ch;counts_fix_ids_ch }

//Check alignment QC
process checkAlignment {

    publishDir params.outDir

    input:
        file(inputBam) from aligned_ch2

    output:
        file("*.RNA_Metrics") into alignCheck_ch

    script:
        """
        java -jar /uufs/chpc.utah.edu/common/home/u0854535/software/picard-2.25.6/picard/build/libs/picard.jar CollectRnaSeqMetrics \
            I=$inputBam \
            O=${inputBam.simpleName}.RNA_Metrics \
            REF_FLAT=$params.refFlat \
            STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
        """ 
        //Do we want to include ribsomal intervals in this check?
        
}


process indexStats {

    publishDir params.outDir

    input:
        path(inputBam) from aligned_ch3
    
    output:
        file("*.stats") into idx_ch
        val true into done_ch

    script:
        """
        ml samtools
        samtools idxstats $inputBam > ${inputBam.simpleName}.stats
        """
}
//R/4.0.2 has DESeq2 installed in its site library (.Library.site) 


process multiqc {

    publishDir params.multi_qc

    input:
        val flag from done_ch

    output:
        file("*.html") into mult_ch

    //Running multiqc from conda working env
    script:
        """
        multiqc $params.outDir $params.bams
        """

}


process deseq2 {

    publishDir params.DESeq

    input:
       path(counts) from counts_ch 

    output:
        file("*") into exp_ch

    script:
        """
        cut -f1,7- $counts > formattedCounts.txt
        ml R/3.6.1
        Rscript $params.projectDir/runDESeq2.R -c formattedCounts.txt -p $params.phenoFile
        """

}


/*
process convertEnsemblIDs {

    publishDir params.outDir

    input:
        path(counts) from counts_fix_ids_ch

    output:
        
    script:


}
*/ 
