#!/usr/bin/env Rscript
# Script for benchmarking import functions for methrix, bsseq, RnBeads, and methylKit

# MIT License
# Copyright (c) 2020 Anand Mayakonda


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-d", default="./bismark_cov/cov_files/",
                    help="Directory containing bismark coverage file [default %(default)s]")
parser$add_argument("-n", type="integer", default=5,
                    help="Number of files to process [default %(default)s]",
                    metavar="number")
parser$add_argument("-t", default="methrix",
                    help = "Tool to use [default \"%(default)s\"]. Can be methrix, methylKit, bsseq, RnBeads")
parser$add_argument("-o", default="./results/import/",
                    help = "Output directory to store results [default \"%(default)s\"]")
parser$add_argument("-f", default=FALSE, action = "store_true",
                    help = "Whether to use ondisk/h5 wherever possible [default \"%(default)s\"].")
parser$add_argument("-c", default=FALSE, action = "store_true",
                    help = "Whether to use CpG loci while running bsseq [default \"%(default)s\"].")


args <- parser$parse_args()

cov_dir = args$d
nfiles = args$n
tool = match.arg(arg = args$t, choices = c("methrix", "bsseq", "methylKit", "RnBeads"))
op_dir= normalizePath(args$o)
h5= args$f
use_cpgs = args$c

#---------------------------------------------------------------------------------------------------------
import_bsseq = function(anno, h5 = FALSE, use_loci = FALSE){
  
  suppressMessages(suppressPackageStartupMessages((library("BiocParallel"))))
  suppressMessages(suppressPackageStartupMessages((library("bsseq"))))
  suppressMessages(suppressPackageStartupMessages((library("data.table"))))
  suppressMessages(suppressPackageStartupMessages((library("BSgenome.Hsapiens.UCSC.hg38"))))
  
  if(use_loci){
    cpg_loci = bsseq::findLoci(pattern = "CG", subject = BSgenome.Hsapiens.UCSC.hg38, strand = "*")  
  }else{
    cpg_loci = NULL
  }
  
  
  if(h5){
    backend = "HDF5Array"
    h5dir = "./bsseq_h5_dir/"
  }else{
    backend = NULL
    h5dir = NULL
  }
  data.table::setDTthreads(threads = 4L)
  res = microbenchmark::microbenchmark(bsseq_obj <- bsseq::read.bismark(anno$bedgraph_files,
                                   loci = cpg_loci,
                                   colData = anno,
                                   rmZeroCov = FALSE,
                                   strandCollapse = TRUE,
                                   BACKEND = backend,
                                   dir = h5dir,
                                   replace = TRUE,
                                   chunkdim = NULL,
                                   level = NULL,
                                   nThread = 4L,
                                   BPPARAM = MulticoreParam(workers = 1L)), times = 1L)
  res
}

#---------------------------------------------------------------------------------------------------------
import_methrix = function(anno, h5 = FALSE){
  
  if(h5){
    backend = TRUE
    h5dir = "./methrix_h5_dir/"
  }else{
    backend = FALSE
    h5dir = NULL
  }
  
  suppressMessages(suppressPackageStartupMessages((library("methrix"))))
  suppressMessages(suppressPackageStartupMessages((library("microbenchmark"))))
  suppressMessages(suppressPackageStartupMessages((library("BSgenome.Hsapiens.UCSC.hg38"))))
  
  hg38_cpgs = suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))
  
  data.table::setDTthreads(threads = 4L)
  res = microbenchmark::microbenchmark(meth_obj = methrix::read_bedgraphs(
    files = anno$bedgraph_files,
    pipeline = "Bismark",
    ref_cpgs = hg38_cpgs,
    zero_based = FALSE,
    stranded = TRUE,
    collapse_strands = TRUE, coldata = sample_anno, h5 = backend, h5_dir = h5dir), times = 1L)
  
  res
}

#---------------------------------------------------------------------------------------------------------
run_methykit_import = function(anno){
  methylKit_obj = methylKit::methRead(
    location = as.list(anno$bedgraph_files),
    sample.id = as.list(anno$Sample_Names),
    assembly = "hg38",
    pipeline = "bismarkCoverage",
    treatment = anno$treatment)
  methylKit_obj <- methylKit::unite(methylKit_obj, destrand=FALSE, min.per.group = 1L)
  methylKit_obj
}

import_methylkit = function(anno){
  suppressMessages(suppressPackageStartupMessages((library("methylKit"))))
  suppressMessages(suppressPackageStartupMessages((library("microbenchmark"))))
  res = microbenchmark::microbenchmark(mobj <- run_methykit_import(anno = anno), times = 1L)
  res
}

#---------------------------------------------------------------------------------------------------------
run_rnbeads = function(anno, dir = "./", h5 = FALSE){
  suppressMessages(suppressPackageStartupMessages(library("RnBeads")))
  suppressMessages(suppressPackageStartupMessages((library("RnBeads.hg38"))))
  suppressMessages(suppressPackageStartupMessages(library("microbenchmark")))
  
  
  anno$bedgraph_files = basename(anno$bedgraph_files)
  RnBeads::rnb.options(import.bed.style = "bismarkCov",
                       region.types = c(),
                       disk.dump.big.matrices = h5, assembly = "hg38")
  dir = normalizePath(dir)
  write.csv(anno, paste0(dir, "/rnbeads.csv"), row.names = FALSE)
  
  setwd(dir)
  
  res = microbenchmark::microbenchmark(rnb.set <-
    rnb.execute.import(data.source = list(dir,
                                          paste0(dir, "/rnbeads.csv"),
                                          2),
                       data.type = "bs.bed.dir"), times = 1L)
  
  res
}
#---------------------------------------------------------------------------------------------------------

bedGraph_files = list.files(cov_dir, pattern = "*\\.bismark\\.cov\\.gz$", full.names = TRUE)
sample_anno = data.frame(
  row.names = gsub(
    pattern = ".bismark.cov.gz",
    replacement = "",
    x = basename(path = bedGraph_files)
  ),
  Sample_Names = gsub(
    pattern = ".bismark.cov.gz",
    replacement = "",
    x = basename(path = bedGraph_files)
  )
)
sample_anno$bedgraph_files = bedGraph_files
sample_anno$Condition = substr(x = sample_anno$Sample_Names, start = 2, stop = 2)
sample_anno$treatment = ifelse(test = sample_anno$Condition == "T", yes = 1, no = 0)

sample_anno = sample_anno[1:nfiles, ]

message("#------------------------------------")
message("Package: ", tool, "\nfiles: ", nfiles)
message("#------------------------------------")

if(tool == "methrix"){
  res = import_methrix(anno = sample_anno, h5 = h5)
  saveRDS(object = res, file = paste0(op_dir,"/", tool, "_", nfiles, "samples", "_", h5, ".RDs"), version = 2)
}else if(tool == "bsseq"){
  res = import_bsseq(anno = sample_anno, h5 = h5, use_loci = use_cpgs)
  if(use_cpgs){
    saveRDS(object = res, file = paste0(op_dir,"/", tool, "Loci_", nfiles, "samples", "_", h5, ".RDs"), version = 2)
  }else{
    saveRDS(object = res, file = paste0(op_dir,"/", tool, "_", nfiles, "samples", "_", h5, ".RDs"), version = 2)   
  }
}else if(tool == "methylKit"){
  res = import_methylkit(anno = sample_anno)
  saveRDS(object = res, file = paste0(op_dir,"/", tool, "_", nfiles, "samples", "_", h5, ".RDs"), version = 2)
}else if(tool == "RnBeads"){
  res = run_rnbeads(anno = sample_anno, dir = cov_dir, h5 = h5)
  saveRDS(object = res, file = paste0(op_dir,"/", tool, "_", nfiles, "samples", "_", h5, ".RDs"), version = 2)
}

message("#------------------------------------")
