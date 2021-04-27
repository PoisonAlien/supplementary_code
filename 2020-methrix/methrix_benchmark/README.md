Scripts used for speed-benchmarking `methrix` `bsseq` `RnBeads` and `methylKit` Bioconductor packages. 

* `runner.sh` contains bash script which runs all benchmarks
* `src` contain individual rscripts for benchmarking spcecific task
* `results` contain output from individual task


### Bismark coverage file sizes

20 Bismark generated coverage files were used for benchmarking. Files were generated as a part of internal project at high coverage using distinct library prparation protocols.


| File                    | n_CpGs   | Library |
|-------------------------|----------|---------|
| 5N_PBAT.bismark.cov.gz  | 45816481 | PBAT    |
| 6N_PBAT.bismark.cov.gz  | 45367792 | PBAT    |
| 5T_PBAT.bismark.cov.gz  | 47224004 | PBAT    |
| 5N_TWGBS.bismark.cov.gz | 53993430 | TWGBS   |
| 5T_TWGBS.bismark.cov.gz | 54975496 | TWGBS   |
| 6T_EMseq.bismark.cov.gz | 54926428 | EMseq   |
| 5N_SWIFT.bismark.cov.gz | 56146529 | SWIFT   |
| 5T_SWIFT.bismark.cov.gz | 54928426 | SWIFT   |
| 5N_EMseq.bismark.cov.gz | 56161262 | EMseq   |
| 6N_TWGBS.bismark.cov.gz | 54422582 | TWGBS   |
| 6N_SWIFT.bismark.cov.gz | 55185050 | SWIFT   |
| 5T_EMseq.bismark.cov.gz | 56259475 | EMseq   |
| 6N_EMseq.bismark.cov.gz | 56288935 | EMseq   |
| 5T_WGBS.bismark.cov.gz  | 56381500 | WGBS    |
| 5N_WGBS.bismark.cov.gz  | 56331423 | WGBS    |
| 6N_WGBS.bismark.cov.gz  | 56454067 | WGBS    |
| 6T_PBAT.bismark.cov.gz  | 46385430 | PBAT    |
| 6T_SWIFT.bismark.cov.gz | 56049895 | SWIFT   |
| 6T_TWGBS.bismark.cov.gz | 54347599 | TWGBS   |
| 6T_WGBS.bismark.cov.gz  | 56517579 | WGBS    |



### SessionInfo

```r
> sessionInfo()
R version 4.0.0 (2020-04-24)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /ngs_share/tools/miniconda3/envs/methrix_benchmark/lib/libopenblasp-r0.3.9.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] methylKit_1.14.1                       
 [2] bsseq_1.24.0                           
 [3] RnBeads_2.6.0                          
 [4] plyr_1.8.6                             
 [5] methylumi_2.34.0                       
 [6] minfi_1.34.0                           
 [7] bumphunter_1.30.0                      
 [8] locfit_1.5-9.4                         
 [9] iterators_1.0.12                       
[10] foreach_1.5.0                          
[11] Biostrings_2.56.0                      
[12] XVector_0.28.0                         
[13] FDb.InfiniumMethylation.hg19_2.2.0     
[14] org.Hs.eg.db_3.11.1                    
[15] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[16] GenomicFeatures_1.40.0                 
[17] AnnotationDbi_1.50.0                   
[18] reshape2_1.4.4                         
[19] scales_1.1.0                           
[20] illuminaio_0.30.0                      
[21] limma_3.44.1                           
[22] gridExtra_2.3                          
[23] gplots_3.0.3                           
[24] ggplot2_3.3.0                          
[25] fields_10.3                            
[26] maps_3.3.0                             
[27] spam_2.5-1                             
[28] dotCall64_1.0-0                        
[29] ff_2.2-14.2                            
[30] bit_1.1-15.2                           
[31] cluster_2.1.0                          
[32] MASS_7.3-51.6                          
[33] methrix_1.2.05                         
[34] SummarizedExperiment_1.18.1            
[35] DelayedArray_0.14.0                    
[36] matrixStats_0.56.0                     
[37] Biobase_2.48.0                         
[38] GenomicRanges_1.40.0                   
[39] GenomeInfoDb_1.24.0                    
[40] IRanges_2.22.1                         
[41] S4Vectors_0.26.0                       
[42] BiocGenerics_0.34.0                    
[43] data.table_1.12.8                      

loaded via a namespace (and not attached):
 [1] BiocFileCache_1.12.0      splines_4.0.0            
 [3] BiocParallel_1.22.0       digest_0.6.25            
 [5] gdata_2.18.0              magrittr_1.5             
 [7] memoise_1.1.0             BSgenome_1.56.0          
 [9] readr_1.3.1               annotate_1.66.0          
[11] R.utils_2.9.2             bdsmatrix_1.3-4          
[13] askpass_1.1               siggenes_1.62.0          
[15] prettyunits_1.1.1         colorspace_1.4-1         
[17] blob_1.2.1                rappdirs_0.3.1           
[19] dplyr_0.8.5               crayon_1.3.4             
[21] RCurl_1.98-1.2            genefilter_1.70.0        
[23] GEOquery_2.56.0           survival_3.1-12          
[25] glue_1.4.0                gtable_0.3.0             
[27] zlibbioc_1.34.0           Rhdf5lib_1.10.0          
[29] HDF5Array_1.16.0          mvtnorm_1.1-0            
[31] DBI_1.1.0                 rngtools_1.5             
[33] Rcpp_1.0.4.6              emdbook_1.3.12           
[35] xtable_1.8-4              progress_1.2.2           
[37] mclust_5.4.6              preprocessCore_1.50.0    
[39] httr_1.4.1                RColorBrewer_1.1-2       
[41] ellipsis_0.3.0            pkgconfig_2.0.3          
[43] reshape_0.8.8             XML_3.99-0.3             
[45] R.methodsS3_1.8.0         dbplyr_1.4.3             
[47] tidyselect_1.0.0          rlang_0.4.6              
[49] munsell_0.5.0             tools_4.0.0              
[51] RSQLite_2.2.0             fastseg_1.34.0           
[53] stringr_1.4.0             bit64_0.9-7              
[55] beanplot_1.2              caTools_1.18.0           
[57] scrime_1.3.5              purrr_0.3.4              
[59] nlme_3.1-147              doRNG_1.8.2              
[61] nor1mix_1.3-0             R.oo_1.23.0              
[63] xml2_1.3.2                biomaRt_2.44.0           
[65] compiler_4.0.0            curl_4.3                 
[67] tibble_3.0.1              stringi_1.4.6            
[69] lattice_0.20-41           Matrix_1.2-18            
[71] permute_0.9-5             multtest_2.44.0          
[73] vctrs_0.2.4               pillar_1.4.4             
[75] lifecycle_0.2.0           bitops_1.0-6             
[77] qvalue_2.20.0             rtracklayer_1.48.0       
[79] R6_2.4.1                  KernSmooth_2.23-17       
[81] codetools_0.2-16          gtools_3.8.2             
[83] assertthat_0.2.1          rhdf5_2.32.0             
[85] openssl_1.4.1             rjson_0.2.20             
[87] withr_2.2.0               GenomicAlignments_1.24.0 
[89] Rsamtools_2.4.0           GenomeInfoDbData_1.2.3   
[91] hms_0.5.3                 quadprog_1.5-8           
[93] coda_0.19-3               tidyr_1.0.3              
[95] base64_2.0                DelayedMatrixStats_1.10.0
[97] bbmle_1.0.23.1            numDeriv_2016.8-1.1
```
