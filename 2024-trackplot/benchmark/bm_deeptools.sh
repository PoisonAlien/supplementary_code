#!/usr/bin/env bash

bed="data/U87/BRD4.bed"
outmat="GSM2634761_U87_H3K4Me3_peaks.mat.gz"

computeMatrix reference-point -R ${bed} --binSize 10 \
-S data/U87/BRD4.bw \
data/U87/BRD4_dBET_2h.bw \
data/U87/BRD4_dBET_24h.bw \
--outFileName ${outmat} --referencePoint center -a 1500 -b 1500 --samplesLabel DMSO dBET_24h dBET_2h \
--numberOfProcessors 1

plotProfile --matrixFile ${outmat} --outFileName ${outmat}.profilePlot.png
plotHeatmap --matrixFile ${outmat} --outFileName ${outmat}.heatmap.png
