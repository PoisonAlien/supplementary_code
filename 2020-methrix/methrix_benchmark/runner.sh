#!/usr/bin/env bash
# Main script calling R functions on bsseq, methrix, methylKit, RnBeads

# MIT License
# Copyright (c) 2020 Anand Mayakonda

mkdir -p results

#activate environment if on vm
#conda activate methrix_benchmark

#Generate objects
objects_dir="results/objects/"
mkdir -p ${objects_dir}
echo "-Generating objects: `date`"
echo "-------------------"
echo "--RnBeads"
./src/generate_objects.R -d ./bismark_cov/cov_files/ -n 10 -t RnBeads -o ${objects_dir} 2> ${objects_dir}/rnbeads.log
echo "--methrix"
./src/generate_objects.R -d ./bismark_cov/cov_files/ -n 10 -t methrix -o ${objects_dir} 2> ${objects_dir}/methrix.log
echo "--methylKit"
./src/generate_objects.R -d ./bismark_cov/cov_files/ -n 10 -t methylKit -o ${objects_dir} 2> ${objects_dir}/methylKit.log
echo "--bsseq"
./src/generate_objects.R -d ./bismark_cov/cov_files/ -n 10 -t bsseq -o ${objects_dir} 2> ${objects_dir}/bsseq.log


#Benchmark SNP filtering
echo "-Running SNP filtering: `date`"
echo "----------------------"
snp_dir="results/snp_filter"
mkdir -p ${snp_dir}
echo "--methrix"
./src/snp_filter.R -o ${objects_dir}/methrix_10samples.RDs 2> ${snp_dir}/methrix.log
echo "--RnBeads"
./src/snp_filter.R -o ${objects_dir}/rnbeads_10samples.RDs 2> ${snp_dir}/rnbeads.log
mv SNPfilter_* ${snp_dir}


#Benchmark DMR subsetting
echo "-Running DMR subsetting: `date`"
echo "-----------------------"
dmrs_dir="./results/DMRs/"
mkdir -p ${dmrs_dir}
echo "--Simulating DMRs"
fetchChromSizes hg38 | head -n 24 > ${dmrs_dir}/hg38.chrom.sizes
./src/simulate_DMRs.R -c ${dmrs_dir}/hg38.chrom.sizes -o ${dmrs_dir} > /dev/null 2> /dev/null
echo "--bsseq"
./src/subset_dmrs.R -o ${objects_dir}/bsseq_10samples.RDs -d ${dmrs_dir} 2> ${dmrs_dir}/bsseq.log
echo "--methrix"
./src/subset_dmrs.R -o ${objects_dir}/methrix_10samples.RDs -d ${dmrs_dir} 2> ${dmrs_dir}/methrix.log
echo "--methylKit"
./src/subset_dmrs.R -o ${objects_dir}/methylKit_10samples.RDs -d ${dmrs_dir} 2> ${dmrs_dir}/methylKit.log
mv subsetDMR*.RDs ${dmrs_dir}


#Benchmark masking 
echo "-Running Mask filters: `date`"
echo "---------------------"
mask_dir="results/mask_filter/"
mkdir -p ${mask_dir}
echo "--methylKit"
./src/mask_filter.R -o ${objects_dir}/methylKitRaw_10samples.RDs 2> ${mask_dir}/methylKit.log
echo "--methrix"
./src/mask_filter.R -o ${objects_dir}/methrix_10samples.RDs 2> ${mask_dir}/methrix.log
echo "--RnBeads"
./src/mask_filter.R -o ${objects_dir}/rnbeads_10samples.RDs 2> ${mask_dir}/RnBeads.log
mv mask_*.RDs ${mask_dir}

#Benchmark import functions 
echo "-Running Import functions: `date`"
echo "-------------------------"
import_dir="results/import/"
mkdir -p ${import_dir}

echo "--In memory import:  `date`"
echo "---bsseq"
echo "----5"
./src/read_bdg.R -t bsseq -n 5 -o results/import/ -d ./bismark_cov/cov_files/ 
echo "----10"
./src/read_bdg.R -t bsseq -n 10 -o results/import/ -d ./bismark_cov/cov_files/ 
echo "----15"
./src/read_bdg.R -t bsseq -n 15 -o results/import/ -d ./bismark_cov/cov_files/ 
echo "----20"
./src/read_bdg.R -t bsseq -n 20 -o results/import/ -d ./bismark_cov/cov_files/ 

echo "---bsseq+CpGLoi"
echo "----5"
./src/read_bdg.R -t bsseq -n 5 -o results/import/ -d ./bismark_cov/cov_files/ -c 
echo "----10"
./src/read_bdg.R -t bsseq -n 10 -o results/import/ -d ./bismark_cov/cov_files/ -c 
echo "----15"
./src/read_bdg.R -t bsseq -n 15 -o results/import/ -d ./bismark_cov/cov_files/ -c 
echo "----20"
./src/read_bdg.R -t bsseq -n 20 -o results/import/ -d ./bismark_cov/cov_files/ -c

echo "---methrix"
echo "----5"
./src/read_bdg.R -t methrix -n 5 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----10"
./src/read_bdg.R -t methrix -n 10 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----15"
./src/read_bdg.R -t methrix -n 15 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----20"
./src/read_bdg.R -t methrix -n 20 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null

echo "---methylKit"
echo "----5"
./src/read_bdg.R -t methylKit -n 5 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----10"
./src/read_bdg.R -t methylKit -n 10 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----15"
./src/read_bdg.R -t methylKit -n 15 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----20"
./src/read_bdg.R -t methylKit -n 20 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null

echo "---RnBeads"
echo "----5"
./src/read_bdg.R -t RnBeads -n 5 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----10"
./src/read_bdg.R -t RnBeads -n 10 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----15"
./src/read_bdg.R -t RnBeads -n 15 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null
echo "----20"
./src/read_bdg.R -t RnBeads -n 20 -o results/import/ -d ./bismark_cov/cov_files/ > /dev/null 2> /dev/null

echo "--On disk import: `date`"
echo "---bsseq"
echo "----5"
./src/read_bdg.R -t bsseq -n 5 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----10"
./src/read_bdg.R -t bsseq -n 10 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----15"
./src/read_bdg.R -t bsseq -n 15 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----20"
./src/read_bdg.R -t bsseq -n 20 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null

echo "---methrix"
echo "----5"
./src/read_bdg.R -t methrix -n 5 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----10"
./src/read_bdg.R -t methrix -n 10 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null 
echo "----15"
./src/read_bdg.R -t methrix -n 15 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----20"
./src/read_bdg.R -t methrix -n 20 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null

echo "---RnBeads"
echo "----5"
./src/read_bdg.R -t RnBeads -n 5 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----10"
./src/read_bdg.R -t RnBeads -n 10 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----15"
./src/read_bdg.R -t RnBeads -n 15 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null
echo "----20"
./src/read_bdg.R -t RnBeads -n 20 -o results/import/ -d ./bismark_cov/cov_files/ -f > /dev/null 2> /dev/null

echo "-Summarizing results: `date`"
./src/summarize.R -i results/import/ -m results/mask_filter/ -d results/DMRs/ -s results/snp_filter
