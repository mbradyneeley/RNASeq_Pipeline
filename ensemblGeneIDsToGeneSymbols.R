#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(httr)
require(httr)
require(jsonlite)

###
# Queries ENSEMBL for the gene symbols via their API, great for featureCounts output
# https://www.biotools.fr/human/ensembl_symbol_converter#
###

###
# The path to your .txt list of ENSEMBL IDs
###

path <- getwd()
IDList <- readLines(paste0(path, "/listOfEnsemblIDs.txt"))

if (length(IDList) == 1) {

###
# A single ID to convert - use a GET request
###

    url = "https://biotools.fr/human/ensembl_symbol_converter/?api=1&id=ENSMUSG00000082108"
    r <- GET(url)
    
    output = fromJSON( content(r, "text"), flatten=TRUE)
    df <- as.data.frame(output)
    df

} else if (length(IDList) > 1) {

###
# Multiple IDs to convert - use a POST request
###

    url = "https://biotools.fr/human/ensembl_symbol_converter/"
    ids = IDList
    ids_json <- toJSON(ids)
    
    body <- list(api=1, ids=ids_json)
    r <- POST(url, body = body)
    
    output = fromJSON( content(r, "text"), flatten=TRUE)
    df <- as.data.frame(output)
    write_tsv(df, "newEnsemblGeneSymbols.txt")

} else {

    print("List either has no IDs or is not named listOfEnsemblIDs.txt")

}
