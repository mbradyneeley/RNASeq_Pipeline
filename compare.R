#!/usr/bin/env Rscript

library("readxl")
library("hciR")
library("optparse")
library("dplyr")
library("tidyverse")

option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, help="path to raw count table", metavar="character")
  #make_option(c("-m", "--metadata"), type="character", default=NULL, help="path to metadata table", metavar="character"),
  #make_option(c("-d", "--design"), type="character", default=NULL, help="path to linear model design file", metavar="character"),
  #make_option(c("-x", "--contrasts_matrix"), type="character", default=NULL, help="path to contrasts matrix file", metavar="character"),
  #make_option(c("-r", "--relevel"), type="character", default=NULL, help="path to factor relevel file", metavar="character"),
  #make_option(c("-k", "--contrasts_list"), type="character", default=NULL, help="path to contrasts list file", metavar="character"),
  #make_option(c("-p", "--contrasts_pairs"), type="character", default=NULL, help="path to contrasts pairs file", metavar="character"),
  #make_option(c("-l", "--genelist"), type="character", default=NULL, help="path to gene list file", metavar="character"),
  #make_option(c("-t", "--logFCthreshold"), type="integer", default=0, help="Log 2 Fold Change threshold for DE genes", metavar="character"),
  #make_option(c("-b", "--batchEffect"), default=FALSE, action="store_true", help="Whether to consider batch effects in the DESeq2 analysis", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate and read input
#if (is.null(opt$counts)){
#  print_help(opt_parser)
#  stop("Counts table needs to be provided!")
#} else {
#  path_count_table = opt$counts
#}

#if (is.null(opt$metadata)){
#  print_help(opt_parser)
#  stop("Metadata table needs to be provided!")
#} else {
#  metadata_path = opt$metadata
#}
#
#if (is.null(opt$design)){
#  print_help(opt_parser)
#  stop("Linear model design file needs to be provided!")
#} else {
#  path_design = opt$design
#}
#
#if(!is.null(opt$relevel)){
#  path_relevel = opt$relevel
#}
#
#if(!is.null(opt$contrasts_matrix)){
#  if(!is.null(opt$contrasts_list) & !is.null(opt$contrasts_pairs)) {
#    stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
#  }
#  path_contrasts_matrix = opt$contrasts_matrix
#}
#
#if(!is.null(opt$contrasts_list)){
#  if(!is.null(opt$contrasts_matrix) & !is.null(opt$contrasts_pairs)) {
#    stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
#  }
#  path_contrasts_list = opt$contrasts_list
#}
#
#if(!is.null(opt$contrasts_pairs)){
#  if(!is.null(opt$contrasts_matrix) & !is.null(opt$contrasts_list)) {
#    stop("Provide only one of contrasts_matrix / contrasts_list / contrasts pairs!")
#  }
#  path_contrasts_pairs = opt$contrasts_pairs
#}
#
#if(!is.null(opt$genelist)){
#  requested_genes_path = opt$genelist
#}

samples <- read_excel("RNAseq_submission_72samples_pheno.xlsx")
samples <- samples[order_samples(samples$Sample_ID), ]
samples$ESRD_SL_5 <- factor(samples$ESRD_SL_5, levels = c(1,0))
samples$eGFR <- round(samples$eGFR, 1)
#samples

#Load the combined featureCounts output.

options(width=110)
# Read in the counts file as a command line arg
# Format tge counts matrix like the counts_marcus.txt,
# it includes geneid
path_count_table <- readr::read_delim('counts_marcus.txt', comment = '#')
#path_count_table <- path_count_table[,7:ncol(path_count_table)]

#Remove 16907 features with zero counts and 18008 features with 10 or fewer reads
#in every sample.

counts <- filter_counts(path_count_table, n = 10)

#In addition, remove 3628 genes with 10 or fewer reads in 90% of samples to create
#a final count matrix with 22121 rows

n <- rowSums(as_matrix(counts) > 10)
sum(n <= 7)
counts <- counts[n > 7,]

#Check genes with the highest number of assigned reads.

options(width=110)
library(hciRdata)
# Had to install the above library after installing the remotes package
# and then used this command... remotes::install_github("HuntsmanCancerInstitute/hciRdata")
n <- rowMeans(as_matrix(counts))
inner_join( dplyr::select(human104, 1:4,8),
  tibble(id= names(n), mean_count = n)) %>%
  mutate(description=trimws(substr(description, 1,36)))  %>%
  arrange(desc(mean_count))

#Run DESeq using ESRD_SL_5 in the design formula  and get the regularized log
#transforms (rlog) for sample visualizations.  These values are similar to log2
#normalized counts except the variance in low count genes is reduced.

dds <- deseq_from_tibble(counts, samples, design = ~ ESRD_SL_5)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 297 genes
# -- DESeq argument 'minReplicatesForReplace' = 7
# -- original counts are preserved in counts(dds)
rld <- r_log(dds)

#Plot the first two principal components using the rlog values from the top 500
#variable genes.

plot_pca(rld, "ESRD_SL_5", tooltip=c("Sample_ID", "STUDY_NO", "eGFR"), ntop=500, width=700)

## Model 2, drop outlier

#Remove the X7 outlier from the sample and count tables and re-run DESeq.

s1 <- filter(samples, Sample_ID != "18992X7")
c1 <- counts[, c(1, which(colnames(counts) %in% s1$Sample_ID))]
dds2 <- deseq_from_tibble(c1, s1, design = ~ ESRD_SL_5)
# -- replacing outliers and refitting for 53 genes
rld2 <- r_log(dds2)

#PC2 likely separates male and female samples.

plot_pca(rld2, "ESRD_SL_5", tooltip=c("Sample_ID", "STUDY_NO", "eGFR"), ntop=500, width=700)

#Compare cases (1) vs. control (0) using a 5% false discover rate (FDR).

res <- results_all(dds2, human104)
# Using adjusted p-value < 0.05
# Adding shrunken fold changes to log2FoldChange
# 1 vs. 0: 0 up and 2 down regulated

#Plot fold changes and p-values in a volcano plot.

plot_volcano(res, pvalue=0.8)

#Plot the mean normalized count and fold change in an MA-plot.

plot_ma(res)

#Cluster 177 genes with an FDR < 25% and scale by
#rows, so values are Z-scores and represent the number of standard deviations
#from the mean rlog value.

x <- top_counts(filter(res, padj < 0.25), rld2, filter=FALSE)
nrow(x)
plot_genes(x, "ESRD_SL_5", scale="row", show_rownames=FALSE, fontsize_col=7,
annotation_names_col=FALSE)

### Save results

#Save the DESeq results and normalized counts to a single Excel
#file in DESeq.xlsx and R objects to a binary data file to load into a new session.

write_deseq(res, dds, rld, human104)
save(dds, rld, res, file="dds.rda")
