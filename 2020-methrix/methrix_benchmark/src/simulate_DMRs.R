#!/usr/bin/env Rscript
# Script for generating simulated DMR of given length

# MIT License
# Copyright (c) 2020 Anand Mayakonda


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-f", type="integer",default=500,
                    help="Minimum number of DMRs [default %(default)s]")
parser$add_argument("-t", type="integer", default=5000,
                    help="Maximum number of DMRs [default %(default)s]")
parser$add_argument("-b", type="integer", default=500,
                    help = "Increment by [default \"%(default)s\"]")
parser$add_argument("-m", type="integer", default=500,
                    help = "Average size of the DMR [default \"%(default)s\"]")
parser$add_argument("-s", type="integer", default=250,
                    help = "SD of the DMR [default \"%(default)s\"]")
parser$add_argument("-o", default="./results/DMRs/",
                    help = "Output directory [default \"%(default)s\"]")
parser$add_argument("-c", default="hg38.chrom.sizes",
                    help = "tsv file with chromosome sizes [default \"%(default)s\"]")


args <- parser$parse_args()

min = args$f
max = args$t
by = args$b
avglen = args$m
sd = args$s
chromsizes = args$c
op_dir= args$o

#---------------------------------------------------------------------------------------------------------

num_dmrs = seq(from = min, to = max, by = by)

lapply(num_dmrs, function(d){
  op_bed = paste0(op_dir, "/DMRs_", d, ".bed")
  file.create(op_bed)
  lapply(abs(rnorm(n = d, mean = avglen, sd = sd)), function(dmrlen){
    cmd = paste0("bedtools random -l ", as.integer(dmrlen) ," -n 1 -g ", chromsizes ,">> ", op_bed)
    system(command = cmd, intern = TRUE)
  })
})
