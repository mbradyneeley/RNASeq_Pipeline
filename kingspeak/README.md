# RNASeq_Pipeline (Combine these READMEs into one)
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
