#!/bin/bash -ue
module load fastqc
fastqc -f fastq gut_1.fq > result1
fastqc -f fastq gut_2.fq > result2
