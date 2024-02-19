#!/usr/bin/env bash

hyperfine --export-csv benchmark/dt_vs_tp_hf.csv --warmup 2 -m 3 -M 5 -n deeptools -n trackplot ./benchmark/profilePlot_deeptools.sh ./benchmark/profilePlot_trackplot.R
