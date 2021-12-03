#!/usr/bin/env bash

#SBATCH --time=06:00:00
#SBATCH -o runDESeq2-%j.out
#SBATCH -e runDESeq2-%j.err
#SBATCH --mail-user=brady.neeley@hsc.utah.edu
#SBATCH --mail-type=END
#SBATCH --account=pezzolesi
#SBATCH --partition=notchpeak

Rscript runDESeq2.R -c counts_marcus.txt -p RNAseq_submission_72samples_pheno.xlsx
