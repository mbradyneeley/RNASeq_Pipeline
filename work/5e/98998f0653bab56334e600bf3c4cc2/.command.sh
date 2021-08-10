#!/bin/bash -ue
module load gcc/8.3.0
STAR --genomeDir /uufs/chpc.utah.edu/common/home/u0854535/RNASeq/scripts/star125          --runThreadN 24          --readFilesIn gut_1.fq gut_2.fq          --twopassMode Basic          --outSAMtype BAM SortedByCoordinate          --quantMode TranscriptomeSAM          --outWigType bedGraph          --outWigStrand Unstranded
