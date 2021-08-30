#!/usr/bin/env nextflow

    //params.reads        = "/scratch/general/lustre/u0854535/RNASeqData/SRA/sra/pairedEndReads/*_{1,2}.fastq.gz"
    params.reads        = "/scratch/general/lustre/u0854535/RNASeqData/SRA/sra/pairedEndReads/*.fastq.gz"
    params.outDir       = "/uufs/chpc.utah.edu/common/home/u0854535/workflows/RNASeq/RNASeq_Pipeline/kingspeak/out"
    params.genomeRef    = "/uufs/chpc.utah.edu/common/home/u0854535/RNASeq/scripts/star125"
    params.refGTF       = "/scratch/general/lustre/u0854535/RNASeqRef/Homo_sapiens.GRCh38.102.gtf"
    //Flat file for hg38/GR38
    params.refFlat      = "/scratch/general/lustre/u0854535/RNASeqRef/refFlat_hg38.txt.gz"
    //TODO:Change the scratch directory name each time you change the scratch dir name
    params.fastqcOut    = "/scratch/general/lustre/u0854535/1frost/results/fastqc"
    params.bams         = "/scratch/general/lustre/u0854535/1frost/results/bam"
    params.multi_qc      = "/scratch/general/lustre/u0854535/1frost/results/multi_qc"

log.info """\
==========================================================================


██████  ███    ██  █████  ███████ ███████  ██████        ███    ██ ███████     
██   ██ ████   ██ ██   ██ ██      ██      ██    ██       ████   ██ ██          
██████  ██ ██  ██ ███████ ███████ █████   ██    ██ █████ ██ ██  ██ █████       
██   ██ ██  ██ ██ ██   ██      ██ ██      ██ ▄▄ ██       ██  ██ ██ ██          
██   ██ ██   ████ ██   ██ ███████ ███████  ██████        ██   ████ ██          
                                              ▀▀                               
                                                                               
██████  ██ ██████  ███████ ██      ██ ███    ██ ███████                        
██   ██ ██ ██   ██ ██      ██      ██ ████   ██ ██                             
██████  ██ ██████  █████   ██      ██ ██ ██  ██ █████                          
██      ██ ██      ██      ██      ██ ██  ██ ██ ██      P E Z Z O L E S I      
██      ██ ██      ███████ ███████ ██ ██   ████ ███████ L A B  -  2 0 2 1 


==========================================================================

reads        = $params.reads
outDir       = $params.outDir
genomeRef    = $params.genomeRef
refGTF       = $params.refGTF
refFlat      = $params.refFlat
fastqcOut    = $params.fastqcOut
bams         = $params.bams
multi_qc     = $params.multi_qc

=========================================================================
"""                

//Check for input reads 
if (params.reads) { Channel.fromFilePairs( "${params.reads}", checkIfExists:true )
    .into { reads_in; reads_in2; reads_in3 }
} else { exit 1, 'Input reads not specified!' }

//Channel.fromFilePairs( "${params.reads}", checkIfExists:true )
//    .into { reads_in; reads_in2; reads_in3 }
    
//reads_in.view()


//TODO: Test removal of optical reads with our RNASeq data from the core
process removeOpticalReads {

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

    tag "$sample_id"

    //publishDir params.outDir    

    input:
        tuple val(sample_id), path(fqs) from trimmedFQs1
        //tuple val(sample_id), path(fqs) from reads_in2
        //The tuple looks like [sample_id, [fastq1, fastq2]]

    //output:
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


//Align reads
process alignReads {

    tag "$sample_id"

    publishDir "${params.bams}", mode: 'copy', pattern: '*.bam'
    //publishDir params.bams

    input:
        //tuple val(sample_id), path(fqs) from sorted_reads_ch
        //tuple val(sample_id), path(fqs) from reads_in3
        tuple val(sample_id), path(fqs) from trimmedFQs2

    output:
        path("*.sortedByCoord.out.bam") into aligned_ch 
    
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
             --outFileNamePrefix ${sample_id} \
             --outTmpKeep All \
             --outStd Log \
             --limitBAMsortRAM 33611120692 \
             --readFilesCommand zcat
        """
}

//Set the bams to a new channel for checking alignment
//Channel.fromPath("${params.bams}/*sortedByCoord.out.bam").into{ aligned_ch3; aligned_ch2; aligned_ch }
aligned_ch.into{ aligned_ch3; aligned_ch2; aligned_ch1 }
//TODO: Need indexes of Bams?
//aligned_ch2.view()
//Count uniquely aligned reads overlapping features in the GTF file
process countFeatures {

    publishDir params.outDir

    beforeScript "echo ${alignedReads}"

    input:
        path(bam) from aligned_ch1

    output:
        //this works file '*.counts.txt' into counts_ch
        path '*.counts.txt' into counts_ch

    script:
        //println "$bam"
        """ 
        featureCounts -T 24 -p -s 2  --largestOverlap -a ${params.refGTF} -o ${bam}.counts.txt $bam
        """
}


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
            O=${inputBam}.RNA_Metrics \
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
        samtools idxstats $inputBam > ${inputBam}.stats
        """
}
//R/4.0.2 has DESeq2 installed in its site library (.Library.site) 


process multiqc {

    publishDir params.multi_qc

    input:
        val flag from done_ch

    output:
        file("*.html") into mult_ch

    script:
        """
        ml multiqc
        multiqc /uufs/chpc.utah.edu/common/home/u0854535/workflows/RNASeq/RNASeq_Pipeline/kingspeak/out /scratch/general/lustre/u0854535/1frost/results/bam
        """

}
