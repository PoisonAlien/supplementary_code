#!/usr/bin/env bash

hyperfine --export-csv benchmark/bioC_vs_tp_hf.csv --warmup 2 -m 3 -M 5 -n Gviz -n karyotyper -n trackplot benchmark/bm_Gviz.R benchmark/bm_karyotypr.R benchmark/bm_trackplot_igv.R
