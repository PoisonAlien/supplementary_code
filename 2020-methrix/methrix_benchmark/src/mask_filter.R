#!/usr/bin/env Rscript
# Script used for benchmarking time-taken for coverage filtering

# MIT License
# Copyright (c) 2020 Anand Mayakonda

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-o", 
                    help="Input object. Required. Can be of methrix, methylRawList or RnBiseqSet")
parser$add_argument("-c", default=5, type="integer",
                    help="coverage threshold. Default 5")

args <- parser$parse_args()

object = args$o
cov = args$c
#---------------------------------------------------------------------------------------------------------

run_methrix = function(obj, cov = 5){
  masking <- microbenchmark::microbenchmark(
    filtered_methrix <- suppressMessages(mask_methrix(obj, low_count = cov, high_quantile = NULL)),
    times=1
  )
  masking
}


run_mkit = function(obj, cov = 5){
  masking <- microbenchmark::microbenchmark(
    filtered_myobj <- methylKit::filterByCoverage(obj, lo.count = cov, lo.perc=NULL,
                                                  hi.count = NULL, hi.perc = NULL),
    times=1
  )
  masking
}

run_rnbeads = function(obj, cov){
  rnb.options(filtering.snp="no",
              filtering.low.coverage.masking = TRUE,
              filtering.coverage.threshold = cov,
              filtering.high.coverage.outliers = TRUE)
  
  masking <- microbenchmark::microbenchmark(
    rnb.set2 <-  rnb.execute.low.coverage.masking(obj, covg.threshold = cov),
    times=1)
  masking
}


message("#------------------------------------")


object = readRDS(object)

if(is(object = object, class2 = "methrix")){
  message("Package: methrix")
  res = run_methrix(obj = object, cov = cov)
  saveRDS(object = res, file = "mask_methrix.RDs")
  message("#microbenchmark took: ", res[[2]]/1e9, " seconds")
}else if(is(object = object, class2 = "methylRawList")){
  message("Package: methylKit")
  res = run_mkit(obj = object, cov = cov)
  saveRDS(object = res, file = "mask_methylKit.RDs")
  message("#microbenchmark took: ", res[[2]]/1e9, " seconds")
}else if(is(object = object, class2 = "RnBiseqSet")){
  message("Package: RnBeads")
  res = run_rnbeads(obj = object, cov = cov)
  saveRDS(object = res, file = "mask_RnBeads.RDs")
  message("#microbenchmark took: ", res[[2]]/1e9, " seconds")
}else{
  stop("Object can only be of class methrix, methylRawList or RnBiseqSet")
  message("#------------------------------------")
}
message("#------------------------------------")
