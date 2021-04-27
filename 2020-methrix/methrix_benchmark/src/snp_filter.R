#!/usr/bin/env Rscript
# Script for benchmarking time-taken by methrix and RnBeads for removing SNPs from objects

# MIT License
# Copyright (c) 2020 Anand Mayakonda

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-o",
                    help="Input object. Can be of methrix, RnBiseqSet class")
parser$add_argument("-m", default=5,type="integer",
                    help = "Minor Allelic Frequency in percent (integer). Default 5 (which will be converted to 0.05).")


args <- parser$parse_args()

object = args$o
maf = args$m/100
#---------------------------------------------------------------------------------------------------------

library(microbenchmark)

filter_methrix_snps = function(obj, maf = 0.05){
  suppressMessages(suppressPackageStartupMessages(library("methrix")))
  suppressMessages(suppressPackageStartupMessages(library("MafDb.1Kgenomes.phase3.GRCh38")))
  obj = methrix::subset_methrix(m = obj, contigs = paste0("chr", 1:22))
  
  res = microbenchmark::microbenchmark(filtered_methrix <- suppressMessages(methrix::remove_snps(m = obj)), times = 1L)
  saveRDS(object = res, file = "SNPfilter_Methrix.RDs")
  message("#microbenchmark took: ", res[[2]]/1e9, " seconds")
  message("#------------------------------------")
}

filter_rnbeads_snps = function(obj){
  suppressMessages(suppressPackageStartupMessages(library("RnBeads")))
  
  res = microbenchmark::microbenchmark(filtered_rnb <- rnb.execute.snp.removal(obj, snp = "any"), times = 1L)
  saveRDS(object = res, file = "SNPfilter_RnBeads.RDs")
  message("#microbenchmark took: ", res[[2]]/1e9, " seconds")
  message("#------------------------------------")
}


object = suppressMessages(suppressWarnings(readRDS(file = object)))

message("#------------------------------------")
if(is(object = object, class2 = "methrix")){
  message("Package: methrix")
  filter_methrix_snps(obj = object, maf = m)
}else if(is(object = object, class2 = "RnBiseqSet")){
  message("Package: RnBeads")
  filter_rnbeads_snps(obj = object)
}else{
  stop("Object can only be of class methrix or RnBiseqSet")
  message("#------------------------------------")
}
