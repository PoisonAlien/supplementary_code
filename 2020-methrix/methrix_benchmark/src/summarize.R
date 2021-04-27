#!/usr/bin/env Rscript
# Script for summarizing results from benchmarking functions

# MIT License
# Copyright (c) 2020 Anand Mayakonda


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-i", default="NULL",
                    help="Directory with import results; from read_bdg.R [default %(default)s]")
parser$add_argument("-m", default="NULL",
                    help="Directory with masking results; from mask_filter.R [default %(default)s]")
parser$add_argument("-d", default="NULL",
                    help="Directory with subsetting results; from subset_dmrs.R [default %(default)s]")
parser$add_argument("-s", default="NULL",
                    help="Directory with SNP filtering results; from snp_filter.R [default %(default)s]")


args <- parser$parse_args()

import_dir = args$i
mask_dir = args$m
subset_dir = args$d
snp_dir = args$s

#---------------------------------------------------------------------------------------------------------

library(ggplot2)
library(data.table)

summarize_import = function(dir = NULL){
  files = list.files(path = dir, pattern = "*\\.RDs$", full.names = TRUE)
  
  message("#----- Summarizing import results")
  message("#----- dir: ", dir)
  if(length(files) == 0 ){
    stop("No files found in ", dir)
  }
  message("#----- nFiles: ", length(files))
  
  
  results = data.frame(data.table::tstrsplit(basename(path = files), split = "_"))
  colnames(results) = c("tool", "n_files", "is_h5")
  results$n_files = paste0("N", gsub(pattern = "samples", replacement = "", x = results$n_files))
  results$tool = ifelse(test = results$is_h5 == "TRUE.RDs", yes = paste0(results$tool, "_h5"), no = results$tool)
  results$file = files
  results$time_secs = unlist(lapply(results$file, function(x) {readRDS(x)[[2]]/1e9}))
  results$tool = factor(x = results$tool, levels = c("bsseq", "bsseqLoci", "bsseq_h5", "methrix", "methrix_h5", "methylKit", "RnBeads", "RnBeads_h5"))
  results$n_files = factor(x = results$n_files, levels = c("N5", "N10", "N15", "N20"))
  results$time_secs_log10 = log10(results$time_secs)
  
  p = ggplot(data = results,
         aes(
           y = time_secs_log10,
           x = n_files,
           group = tool,
           color = tool,
           label = tool
         )) + geom_line() + scale_colour_brewer(palette = "Dark2") + geom_point() +
    theme(text = element_text(size = 16))
  
  pdf(file = paste0(dir, "/import_time_taken.pdf"), width = 8, height = 5, bg = "white", paper = "special")
  print(p)
  dev.off()
  
  p2 = ggplot(data = results,
             aes(
               y = time_secs,
               x = n_files,
               group = tool,
               color = tool,
               label = tool
             )) + geom_line() + scale_colour_brewer(palette = "Dark2") + geom_point() +
    theme(text = element_text(size = 16))
  
  pdf(file = paste0(dir, "/import_time_taken_log10.pdf"), width = 8, height = 5, bg = "white", paper = "special")
  print(p2)
  dev.off()
  
  data.table::fwrite(x = results,  paste0(dir, "/import_time_taken.tsv"), sep = "\t")
  message("#----- Done!")
}


summarize_snp = function(dir = NULL){
  files = list.files(path = dir, pattern = "SNPfilter_*.*\\.RDs$", full.names = TRUE)
  
  if(length(files) == 0 ){
    stop("No files found in ", dir)
  }
  
  message("#----- Summarizing SNP results")
  message("#----- dir: ", dir)
  if(length(files) == 0 ){
    stop("No files found in ", dir)
  }
  message("#----- nFiles: ", length(files))
  
  results = data.frame(data.table::tstrsplit(basename(path = files), split = "_", keep = 2))
  colnames(results) = c("tool")
  results$file = files
  results$tool = gsub(pattern = "\\.RDs", replacement = "", x = results$tool)
  results$time_secs = unlist(lapply(results$file, function(x) {readRDS(x)[[2]]/1e9}))
  
  pdf(file = paste0(dir, "/SNP_filtering.pdf"), width = 4, height = 5, bg = "white", paper = "special")
  barplot(results$time_secs, names.arg = results$tool, ylab = "Time (secs)", ylim = c(0, max(results$time_secs)))
  title(main = "SNP filtering", sub = c("methrix::remove_snps()\nRnBeads::rnb.execute.snp.removal()"), adj = 0, font.sub = 3, cex.sub = 1)
  dev.off()
  data.table::fwrite(x = results,  paste0(dir, "/SNPfiltering_time_taken.tsv"), sep = "\t")
  message("#----- Done!")
}

summarize_dmrs = function(dir = NULL){
  files = list.files(path = dir, pattern = "subsetDMR_*.*\\.RDs$", full.names = TRUE)
  
  if(length(files) == 0 ){
    stop("No files found in ", dir)
  }
  
  message("#----- Summarizing subsetting results")
  message("#----- dir: ", dir)
  if(length(files) == 0 ){
    stop("No files found in ", dir)
  }
  message("#----- nFiles: ", length(files))
  
  results = data.frame(data.table::tstrsplit(basename(path = files), split = "_", keep = 2))
  colnames(results) = c("tool")
  results$file = files
  results$tool = gsub(pattern = "\\.RDs", replacement = "", x = results$tool)
  results_dat = lapply(results$file, function(x) {
    xdat = readRDS(x)
    xdat = lapply(xdat, function(y) y[[2]]/1e9)
    xdat = data.frame(time_secs = unlist(xdat))
    xdat$n_files = paste0("N",
           gsub(
             pattern = ".bed",
             replacement = "",
             x = data.table::tstrsplit(rownames(xdat), split = "_")[[2]]
           ))
    data.table::setDT(xdat, keep.rownames = TRUE)
    xdat
    })
  names(results_dat) = results$tool
  results_dat = data.table::rbindlist(l = results_dat, idcol = "tool")
  results_dat$n_files = factor(x = results_dat$n_files, levels = paste0("N", seq(500, 5000, 500)))
  
  p = ggplot(data = results_dat, aes(x = tool, y = time_secs)) + geom_bar(stat = "identity") +
    facet_grid( ~ n_files) + theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )) + ggtitle(label = "Subset", subtitle = "bsseq::subsetByOverlaps()\nmethrix::subset_methrix\nmethylKit::selectByOverlap()")
  
  
  pdf(file = paste0(dir, "/subsetting_time_taken.pdf"), width = 7, height = 6, bg = "white", paper = "special")
  print(p)
  dev.off()
  
  data.table::fwrite(x = results_dat,  paste0(dir, "/subsetting_time_taken.tsv"), sep = "\t")
  message("#----- Done!")
}


summarize_mask = function(dir = NULL){
  files = list.files(path = dir, pattern = "mask_*.*\\.RDs$", full.names = TRUE)
  
  message("#----- Summarizing masking results")
  message("#----- dir: ", dir)
  if(length(files) == 0 ){
    stop("No files found in ", dir)
  }
  message("#----- nFiles: ", length(files))
  
  results = data.frame(data.table::tstrsplit(basename(path = files), split = "_", keep = 2))
  colnames(results) = c("tool")
  results$file = files
  results$tool = gsub(pattern = "\\.RDs", replacement = "", x = results$tool)
  results$time_secs = unlist(lapply(results$file, function(x) {readRDS(x)[[2]]/1e9}))
  
  pdf(file = paste0(dir, "/Masking_time_taken.pdf"), width = 5, height = 6, bg = "white", paper = "special")
  barplot(results$time_secs, names.arg = results$tool, ylab = "Time (secs)", ylim = c(0, max(results$time_secs)))
  title(main = "Masking low-coverage CpGs", sub = c("methrix::mask_methrix()\nmethylKit::filterByCoverage()\nRnBeads::rnb.execute.low.coverage.masking()"), adj = 0, font.sub = 3, cex.sub = 1)
  dev.off()
  data.table::fwrite(x = results,  paste0(dir, "/Masking_time_taken.tsv"), sep = "\t")
  message("#----- Done!")
}

#---------------------------------------------------------------------------------------------------------

if(import_dir != "NULL"){
  summarize_import(dir = import_dir)  
}

if(snp_dir != "NULL"){
  summarize_snp(dir = snp_dir)
}

if(subset_dir != "NULL"){
  summarize_dmrs(dir = subset_dir)
}

if(mask_dir != "NULL"){
  summarize_mask(dir = mask_dir)
}

