#!/usr/bin/env Rscript

library(data.table)
suppressMessages(suppressPackageStartupMessages((library("methrix"))))
suppressMessages(suppressPackageStartupMessages((library("BSgenome.Hsapiens.UCSC.hg38"))))

hg38_cpgs = suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))

#---------------------------------------------------------------------------------------------------------

import_methrix = function(anno, cpgs = NULL, nbatch = 1, threads = 2){
  
  data.table::setDTthreads(threads = 4L)
  message("Methrix: Importing..")
  meth_obj = methrix::read_bedgraphs(
    files = anno$bedgraph_files,
    pipeline = "Bismark",
    ref_cpgs = cpgs,
    zero_based = FALSE,
    stranded = TRUE,
    collapse_strands = TRUE, coldata = anno, h5 = FALSE, vect = TRUE, n_threads = threads)
  
  meth_obj
}

#---------------------------------------------------------------------------------------------------------

bedGraph_files = list.files(path = "bismark_cov/", pattern = ".cov.gz", recursive = TRUE, full.names = TRUE)
sample_anno = data.frame(source = basename(dirname(bedGraph_files)))
sample_anno$bedgraph_files = bedGraph_files
sample_anno$Condition = substr(x = gsub(pattern = ".bismark.cov.gz", replacement = "", x = basename(sample_anno$bedgraph_files)), start = 2, stop = 2)
sample_anno$treatment = ifelse(test = sample_anno$Condition == "T", yes = 1, no = 0)
rownames(sample_anno) = paste(sample_anno$source, gsub(pattern = ".bismark.cov.gz", replacement = "", x = basename(sample_anno$bedgraph_files)), sep = "_")

#---------------------------------------------------------------------------------------------------------

methrix_parallel_bm = microbenchmark::microbenchmark(
  import_methrix(anno = sample_anno, cpgs = hg38_cpgs, threads = 2),
  import_methrix(anno = sample_anno, cpgs = hg38_cpgs, threads = 4),
  import_methrix(anno = sample_anno, cpgs = hg38_cpgs, threads = 6),
  import_methrix(anno = sample_anno, cpgs = hg38_cpgs, threads = 8),
  times = 3
)

saveRDS(object = methrix_parallel_bm, file = "01_bedgraph_import/methrix_parallel_bm_nthreads.RDs", version = 2)

methrix_parallel_bm$n_threads = unlist(lapply(as.character(methrix_parallel_bm$expr), function(x) {xthr = unlist(data.table::tstrsplit(x = x, split = ",", keep = 3)); xthr = substr(x = xthr, start = 2, stop = nchar(xthr)-1)}))
prll_summary = data.frame(threads = methrix_parallel_bm$n_threads, time = methrix_parallel_bm$time/1e9/60)
prll_summary$threads = as.factor(x = bm_summary$threads)

pdf(file = "bm_summary.pdf", width = 3, height = 4, paper = "special", bg = "white")
par(mar = c(3, 3, 1, 1))
b = boxplot(time ~ threads, data = prll_summary, axes = FALSE, ylim = c(0, 8), xlab = NA, ylab = NA)
axis(side = 2, at = pretty(seq(0, 8)), las = 2, las = 2)
axis(side = 1, at = 1:4, labels = seq(2, 8, 2))
text(x = 1:4, y = apply(b$stats, 2, max)+0.25, labels = round(apply(b$stats, 2, max), digits = 2))
mtext(text = "cores (N)", side = 1, line = 2)
mtext(text = "time (min)", side = 2, line = 2)
dev.off()