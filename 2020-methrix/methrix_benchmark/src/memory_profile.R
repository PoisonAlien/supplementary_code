setwd("/home/anand/Documents/methrix_benchmark")
bedGraph_files = list.files("bismark_cov/cov_files/",
                            pattern = ".bismark.cov.gz",
                            full.names = TRUE)
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
sample_anno$Condition = substr(x = sample_anno$Sample_Names,
                               start = 2,
                               stop = 2)
sample_anno$treatment = ifelse(test = sample_anno$Condition == "T",
                               yes = 1,
                               no = 0)

sample_anno = sample_anno[1:2, , drop = FALSE]

#--------------------------------------------------------------------------------------------
library("BSgenome.Hsapiens.UCSC.hg38")
library("BiocParallel")
cpgs = bsseq::findLoci(pattern = "CG",
                       subject = BSgenome.Hsapiens.UCSC.hg38,
                       strand = "*")
library("bsseq")

bsseq_profvis = profvis::profvis(expr =  {
  bsseq_obj <- bsseq::read.bismark(
    sample_anno$bedgraph_files,
    loci = cpgs,
    colData = sample_anno,
    rmZeroCov = FALSE,
    strandCollapse = TRUE,
    BACKEND = NULL,
    dir = h5dir,
    replace = TRUE,
    chunkdim = NULL,
    level = NULL,
    nThread = 4L,
    BPPARAM = MulticoreParam(workers = 1L)
  )
}, prof_output = "bsseq_read")

saveRDS(object = bsseq_profvis, file = "bsseq.RProf")
htmlwidgets::saveWidget(bsseq_profvis, "bsseq.RProf.html")
gc()

#--------------------------------------------------------------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")
library("methrix")
hg38_cpgs = suppressWarnings(methrix::extract_CPGs(ref_genome = "BSgenome.Hsapiens.UCSC.hg38"))

methrix_profvis = profvis::profvis(expr = {meth_obj <- methrix::read_bedgraphs(
  files = sample_anno$bedgraph_files,
  pipeline = "Bismark",
  ref_cpgs = hg38_cpgs,
  zero_based = FALSE,
  stranded = TRUE,
  collapse_strands = TRUE,
  coldata = sample_anno,
  h5 = FALSE)})

saveRDS(object = methrix_profvis, file = "methrix.RProf.RDs")
htmlwidgets::saveWidget(methrix_profvis, "methrix.RProf.html")
gc()

#--------------------------------------------------------------------------------------------

library("methylKit")
run_mkit = function(anno){
  methylKit_obj = methylKit::methRead(
    location = as.list(anno$bedgraph_files),
    sample.id = as.list(anno$Sample_Names),
    assembly = "hg38",
    pipeline = "bismarkCoverage",
    treatment = anno$treatment
  )
  methylKit_obj <-
    methylKit::unite(methylKit_obj, destrand = FALSE, min.per.group = 1L)
  methylKit_obj
}

methylkit_profvis = profvis::profvis(expr = {
  run_mkit(anno = sample_anno)
}, prof_output = "methylkit.RProf")

saveRDS(object = methylkit_profvis, file = "methylkit.RProf.RDs")
htmlwidgets::saveWidget(methylkit_profvis, "methylkit.RProf.html")
gc()

#--------------------------------------------------------------------------------------------
suppressMessages(suppressPackageStartupMessages(library("RnBeads")))
suppressMessages(suppressPackageStartupMessages((library("RnBeads.hg38"))))

run_rnbeads = function(anno, h5 = FALSE){
  anno$bedgraph_files = basename(anno$bedgraph_files)
  RnBeads::rnb.options(import.bed.style = "bismarkCov",
                       region.types = c(),
                       disk.dump.big.matrices = h5, assembly = "hg38")

  write.csv(anno, "./bismark_cov/cov_files/rnbeads.csv", row.names = FALSE)

  setwd("./bismark_cov/cov_files/")

  rnb.set <- rnb.execute.import(
    data.source = list(
      "/home/anand/Documents/methrix_benchmark/bismark_cov/cov_files/",
      "/home/anand/Documents/methrix_benchmark/bismark_cov/cov_files/rnbeads.csv",
      2
    ),
    data.type = "bs.bed.dir"
  )
  rnb.set
}

rnbeads_profvis = profvis::profvis(expr = {
  run_rnbeads(anno = sample_anno)
}, prof_output = "RnBeads.RProf")

saveRDS(object = rnbeads_profvis, file = "rnbeads.RProf.RDs")
htmlwidgets::saveWidget(rnbeads_profvis, "rnbeads.RProf.html")
gc()
