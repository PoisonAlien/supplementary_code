#!/usr/bin/env Rscript
# Script used for benchmarking time-taken for coverage filtering

# MIT License
# Copyright (c) 2020 Anand Mayakonda

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-b", 
                    help="Input bsseq object. Required.")
parser$add_argument("-r", 
                    help="Input RnBeads object. Required.")
parser$add_argument("-m", 
                    help="Input methrix object. Required.")
parser$add_argument("-c", default=5,type="integer",
                    help="Coverage threshold. Default 5.")
parser$add_argument("-s", default=2,type="integer",
                    help="Minimum number of sampels with coverage. Default 2.")



args <- parser$parse_args()
#---------------------------------------------------------------------------------------------------------

library(bsseq)
library(methrix)
library(RnBeads)
library(microbenchmark)

methrix_obj = readRDS(args$m)
if(is(object = methrix_obj, class2 = "methrix")){
  stop(args$r, "is not a methrix object!") 
}

bsseq_obj = readRDS(args$b)
if(is(object = bsseq_obj, class2 = "bsseq")){
  stop(args$m, "is not a bsseq object!") 
}

rnb_obj = readRDS(args$r)
if(is(object = rnb_obj, class2 = "RnBiseqSet")){
  stop(args$r, "is not a RnBiseqSet object!") 
}

cov_filt_bsseq <- function(bsseq_obj, cov_thr = 5, min_samples = 2){
  BS.cov <- bsseq::getCoverage(bsseq_obj)
  keepLoci.ex <- which(rowSums(BS.cov >= cov_thr) > min_samples)
  bsseq_obj_cov <- bsseq_obj[keepLoci.ex,]
}

cov_filt <- microbenchmark::microbenchmark(
  filtered_methrix <- suppressMessages(methrix::coverage_filter(m = methirx_obj, cov_thr = args$c, min_samples = args$s)),
  filtered_bsseq <- cov_filt_bsseq(bsseq_obj = bsseq_obj, cov_thr = args$c, min_samples = args$s),
  times=1
)
