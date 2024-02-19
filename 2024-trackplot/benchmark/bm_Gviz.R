#!/usr/bin/env Rscript

#load lib

library(Gviz)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#ideogram
#itrack <- Gviz::IdeogramTrack(genome = "hg19", chromosome = "chr6", name = "chr6")


from <- 31125776
to <- 31144789
chr = "chr6"

#gene model
#gtrack <- GeneRegionTrack(chromosome = chr, start = from, end = to, genome = "hg19", range = TxDb.Hsapiens.UCSC.hg19.knownGene)

#Axis track
atrack <- Gviz::GenomeAxisTrack()

#data track
h1_bigWigs = list.files(path = "/Users/anandmayakonda/Documents/Documents_MacBookPro_work/github/trackplot_bm.nosync/data/H1/", pattern = "\\.bw$", full.names = TRUE)
h1_dtracks = lapply(h1_bigWigs, function(bw){
  Gviz::DataTrack(start = from, end = to, chromosome = "chr6", range = bw, genome = "hg19", name = basename(bw))  
})


#tracklist
tracklist = c(h1_dtracks, list(atrack))

Gviz::plotTracks(trackList = tracklist, from = 31125776, to = 31144789, type = "histogram")
