#!/usr/bin/env Rscript

source("trackplot.R")

bws = c("data/U87/BRD4.bw", "data/U87/BRD4_dBET_2h.bw", "data/U87/BRD4_dBET_24h.bw")
bed = "data/U87/BRD4.bed.gz"

bws = read_coldata(bws = bws, sample_names = c("DMSO", "dBET_24h", "dBET_2h"), build = "hg19")
#Center and extend 1500 both ways from the peak center. Estimate mean signal
pe_bed = profile_extract(colData = bws, bed = bed, startFrom = "center", up = 1500, down = 1500, nthreads = 1)
pe_bed_sum = profile_summarize(sig_list = pe_bed) 

#Profile plot
png("profile_plot.png")
profile_plot(pe_bed_sum)
dev.off()

#Heatmap
png("heatmap_plot.png")
profile_heatmap(pe_bed, zmaxs = 0.8)
dev.off()