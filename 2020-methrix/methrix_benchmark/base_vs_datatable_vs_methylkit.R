#!/usr/bin/env Rscript

library(methylKit)
library(data.table)
suppressMessages(suppressPackageStartupMessages((library("methrix"))))
suppressMessages(suppressPackageStartupMessages((library("BSgenome.Hsapiens.UCSC.hg38"))))

hg38_cpgs = suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))

#---------------------------------------------------------------------------------------------------------

run_methykit_import = function(anno, nthreads = 4){
  message("MethylKit: Importing..")
  methylKit_obj = methylKit::methRead(
    location = as.list(anno$bedgraph_files),
    sample.id = as.list(anno$Sample_Names),
    assembly = "hg38",
    pipeline = "bismarkCoverage",
    treatment = anno$treatment, mincov = 1)
  
  methylKit_obj <- methylKit::unite(methylKit_obj, destrand=FALSE, min.per.group = 1L, mc.cores = nthreads)
  methylKit_obj
}

import_methylkit = function(anno){
  suppressMessages(suppressPackageStartupMessages((library("methylKit"))))
  suppressMessages(suppressPackageStartupMessages((library("microbenchmark"))))
  mb_res = microbenchmark::microbenchmark(mobj <- run_methykit_import(anno = anno), times = 1L)
  mb_res
}

#---------------------------------------------------------------------------------------------------------

run_dt = function(anno, nthreads = 4){
  library(data.table)
  data.table::setDTthreads(threads = nthreads)
  
  message("data.table: Reading files..")
  bdg_dt = lapply(anno$bedgraph_files, function(bdg){
    data.table::fread(file = bdg,
                      nThread = nthreads,
                      sep = "\t",
                      col.names = c("chr", "start", "end", "Beta", "M", "U"),
                      key = c("chr", "start", "end")
    )
  })
  
  names(bdg_dt) = basename(anno$bedgraph_files)
  
  message("data.table: Merging..")
  merged_data = bdg_dt[[1]][,.(chr, start, Beta, M)]
  colnames(merged_data)[3:4] = paste0(c("Beta_", "M_"), names(bdg_dt)[1])
  
  for(idx in 2:length(bdg_dt)){
    merged_data = merge(merged_data, bdg_dt[[idx]][,.(chr, start, Beta, M)], by = c("chr", "start"), all = TRUE)
    colnames(merged_data)[rev(rev(seq_len(ncol(merged_data)))[1:2])] = paste0(c("Beta_", "M_"), names(bdg_dt)[idx])
  }
  rm(bdg_dt)
  
  merged_data
}

#---------------------------------------------------------------------------------------------------------

run_base = function(anno, nthreads = 4){
  
  message("Base: Reading files..")
  bdg_dt = parallel::mclapply(anno$bedgraph_files, function(bdg){
    read.delim(file = bdg, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "Beta", "M", "U"), stringsAsFactors = FALSE,
                   colClasses = c("character", "numeric", "numeric", "numeric", "numeric", 'numeric'), as.is = TRUE)
  }, mc.cores = nthreads)
  names(bdg_dt) = basename(anno$bedgraph_files)
  
  message("Base: Merging..")
  merged_data = bdg_dt[[1]][,c("chr", "start", "Beta", "M")]
  colnames(merged_data)[3:4] = paste0(c("Beta_", "M_"), names(bdg_dt)[1])
  
  for(idx in 2:length(bdg_dt)){
    merged_data = merge(merged_data, bdg_dt[[idx]][,c("chr", "start", "Beta", "M")], by = c("chr", "start"), all = TRUE)
    colnames(merged_data)[rev(rev(seq_len(ncol(merged_data)))[1:2])] = paste0(c("Beta_", "M_"), names(bdg_dt)[idx])
  }
  rm(bdg_dt)
  
  merged_data
}
#---------------------------------------------------------------------------------------------------------

import_methrix = function(anno, cpgs = NULL, nthreads = 4){
  
  data.table::setDTthreads(threads = nthreads)
  message("Methrix: Importing..")
  meth_obj = methrix::read_bedgraphs(
    files = anno$bedgraph_files,
    pipeline = "Bismark",
    ref_cpgs = cpgs,
    zero_based = FALSE,
    stranded = TRUE,
    collapse_strands = TRUE, coldata = anno, h5 = FALSE, n_threads = nthreads)
  
  meth_obj
}

#---------------------------------------------------------------------------------------------------------

bedGraph_files = list.files("bismark_cov/cov_files/", pattern = ".bismark.cov.gz", full.names = TRUE)
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

#Base vs Data.table vs MethylKit
base_dt_mkit_methrix = microbenchmark::microbenchmark(
  run_base(anno = sample_anno, nthreads = 4),
  run_dt(anno = sample_anno, nthreads = 4),
  run_methykit_import(anno = sample_anno, nthreads = 4),
  import_methrix(anno = sample_anno, cpgs = hg38_cpgs, nthreads = 4),
  times = 3
)

saveRDS(object = base_dt_mkit_methrix, file = "01_bedgraph_import/base_dt_mkit_methrix.RDs", version = 2)

base_dt_mkit_methrix$tool = unlist(lapply(as.character(base_dt_mkit_methrix$expr), function(x) unlist(data.table::tstrsplit(x = x, split = "\\(", keep = 1))))
bm_summary = data.frame(tool = base_dt_mkit_methrix$tool, time = base_dt_mkit_methrix$time/1e9/60)
bm_summary = bm_summary[!bm_summary$tool %in% "run_methykit_import",] #Not including methylKit
bm_summary$tool = factor(x = bm_summary$tool, levels = c("run_base", "run_dt", "import_methrix"), labels = c("base", "data.table", "methrix"))

pdf(file = "methrix_DT_mkit_summary.pdf", width = 3, height = 4, paper = "special", bg = "white")
par(mar = c(3, 3, 1, 1))
b = boxplot(time ~ tool, data = bm_summary, axes = FALSE, ylim = c(0, 30), xlab = NA, ylab = NA)
axis(side = 2, at = pretty(seq(0, 30)), las = 2, las = 2)
axis(side = 1, at = 1:3, labels = c("base", "data.table", "methrix"))
text(x = 1:3, y = apply(b$stats, 2, max)+1, labels = round(apply(b$stats, 2, max), digits = 2))
mtext(text = "cores (N)", side = 1, line = 2)
mtext(text = "time (min)", side = 2, line = 2)
dev.off()