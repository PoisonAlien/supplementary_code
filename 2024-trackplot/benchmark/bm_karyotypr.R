#!/usr/bin/env Rscript

library(karyoploteR)

region <- toGRanges("chr6:31125776-31144789")
kp <- karyoploteR::plotKaryotype(zoom = region)

#genes.data <- karyoploteR::makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot=kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)

#kp <- karyoploteR::plotKaryotype(zoom = region)

#genes.data <- addGeneNames(genes.data)
#genes.data <- mergeTranscripts(genes.data)

#kp <- plotKaryotype(zoom = region, cex=1)
#kpPlotGenes(kp, data=genes.data)

#kp <- plotKaryotype(zoom = region, cex=2)
#kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 0.8)

h1_bigWigs = list.files(path = "/Users/anandmayakonda/Documents/Documents_MacBookPro_work/github/trackplot_bm.nosync/data/H1/", pattern = "\\.bw$", full.names = TRUE)
for(i in seq_len(length(h1_bigWigs))) {
  bigwig.file <- h1_bigWigs[i]
  at <- autotrack(i, length(h1_bigWigs), r0=0.35, r1=1, margin = 0.1)
  kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region", 
                     r0=at$r0, r1=at$r1)
  computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1)
  kpAddLabels(kp, labels = names(h1_bigWigs)[i], r0=at$r0, r1=at$r1, 
              cex=1.6, label.margin = 0.035)
}
