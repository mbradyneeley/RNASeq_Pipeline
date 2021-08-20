#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH -o testLonepeak-%j.out
#SBATCH -e testLonepeak-%j.err
#SBATCH --mail-user=brady.neeley@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi-np
#SBATCH --partition=pezzolesi-np
#SBATCH --mem=144G

fq1=/scratch/general/lustre/u0854535/RNASeqData/JC0180_1.fastq.gz
fq2=/scratch/general/lustre/u0854535/RNASeqData/JC0180_2.fastq.gz


module load gcc/8.3.0

STAR --genomeDir /uufs/chpc.utah.edu/common/home/u0854535/RNASeq/scripts/star125 \
     --runThreadN 24 \
     --readFilesIn $fq1 $fq2 \
     --twopassMode Basic \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outWigType bedGraph \
     --outWigStrand Unstranded \
     --outFileNamePrefix ./test/JC0180 \
     --outTmpKeep All \
     --outStd Log \
     --readFilesCommand zcat
