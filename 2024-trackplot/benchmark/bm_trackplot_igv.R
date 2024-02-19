#!/usr/bin/env Rscript

library(trackplot)

h1_bigWigs = list.files(path = "/Users/anandmayakonda/Documents/Documents_MacBookPro_work/github/trackplot_bm.nosync/data/H1/", pattern = "\\.bw$", full.names = TRUE)
h1_bigWigs = trackplot::read_coldata(bws = h1_bigWigs, build = "hg19")

#Region to plot
oct4_loci = "chr6:31125776-31144789"

#Extract bigWig signal for a loci of interest
t = track_extract(colData = h1_bigWigs, loci = oct4_loci, query_ucsc = FALSE)

trackplot::track_plot(
  summary_list = t, groupAutoScale = FALSE, y_min = 0, draw_gene_track = FALSE)
