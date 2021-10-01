# RNASeq_Pipeline
This is a fairly simple pipeline for RNASeq data analysis

Important components of this pipeline are the correct reference made with STAR
genomeGenerate and the refFlat file (solely for QC purposes) used by Picard's 
CollectRnaSeqMetrics. Make your own refFlat file using gtfToGenePred.

The command I used was:

gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    -ignoreGroupsWithoutExons \
    <path/to/your/GTF> \
    /dev/stdout | \
    awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'

and written to a file named refFlat_<GTFBuildVersion>.txt then gzipped.

ToDo:
- Implement value channels in place of some queue channels
    - https://seqera.io/training/
- path() is preferred over file() in the most recent nextflow v.
- Multiqc runs for every one of the input files, should only run once after
    all alignment steps are complete

Requirements:
    -   multiqc

To view logs (especially for lengthy alignment processes) find the current
work directory in theBifrost scratch directory and cat .command.log. Then
go to the given scratch working directory for that process, usually looks
like /scratch/general/lustre/u0854535/temp/nxf.E8161SdX0U
