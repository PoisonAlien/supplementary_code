#!/usr/bin/env Rscript
# Script used for benchmarking time-taken for coverage filtering

# MIT License
# Copyright (c) 2020 Anand Mayakonda

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-o", 
                    help="Input object. Required. Can be of methrix, BSseq or methylBase object")
parser$add_argument("-d", default="./DMRs/",
                    help="Directory containing bed files. Required. All \\.BED$ files are greped")

args <- parser$parse_args()

object = args$o
dir = args$d
#---------------------------------------------------------------------------------------------------------

subset_mkit = function(object, dmr_files){
  suppressMessages(suppressPackageStartupMessages(library("methylKit")))
  tt = lapply(dmr_files, function(x){
    gr = data.table::fread(file = x)
    gr = GenomicRanges::GRanges(seqnames = gr$V1, ranges = IRanges(start = gr$V2, end = gr$V3), strand = gr$V6)
    microbenchmark::microbenchmark(methylKit::selectByOverlap(object = object, gr), times = 1L)
  })
  names(tt) = basename(dmr_files)
  tt
}

subset_methrix = function(object, dmr_files){
  suppressMessages(suppressPackageStartupMessages(library("methrix")))
  tt = lapply(dmr_files, function(x){
    gr = data.table::fread(file = x)
    microbenchmark::microbenchmark(methrix::subset_methrix(m = object, regions =
                                                             gr), times = 1L)
  })
  names(tt) = basename(dmr_files)
  tt
}

subset_bsseq = function(object, dmr_files){
  suppressMessages(suppressPackageStartupMessages(library("bsseq")))
  tt = lapply(dmr_files, function(x){
    gr = data.table::fread(file = x)
    gr = GenomicRanges::GRanges(seqnames = gr$V1, ranges = IRanges(start = gr$V2, end = gr$V3), strand = gr$V6)
    microbenchmark::microbenchmark(bsseq::subsetByOverlaps(object, ranges = gr), times = 1L)
  })
  names(tt) = basename(dmr_files)
  tt
}

bed_files = list.files(dir, pattern = "*\\.bed$", full.names = TRUE)

if(length(bed_files) == 0){
  stop("No BED files found in ", dir)
}

object = suppressPackageStartupMessages(suppressMessages(suppressWarnings(readRDS(file = object))))

message("#------------------------------------")
if(is(object = object, class2 = "methrix")){
  message("Package: methrix", "\nobject: ", args$o, "\nnDMRs: ", length(bed_files))
  res = subset_methrix(obj = object, dmr_files = bed_files)
  saveRDS(object = res, file = "subsetDMR_methrix.RDs")
}else if(is(object = object, class2 = "BSseq")){
  message("Package: bsseq", "\nobject: ", args$o, "\nnDMRs: ", length(bed_files))
  res = subset_bsseq(obj = object, dmr_files = bed_files)
  saveRDS(object = res, file = "subsetDMR_bsseq.RDs")
}else if(is(object = object, class2 = "methylBase")){
  message("Package: methylKit", "\nobject: ", args$o, "\nnDMRs: ", length(bed_files))
  res = subset_mkit(obj = object, dmr_files = bed_files)
  saveRDS(object = res, file = "subsetDMR_methylKit.RDs")
}else{
  stop("Object can only be of class methrix, BSseq or methylBase")
  message("#------------------------------------")
}
message("#------------------------------------")