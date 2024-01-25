#!/usr/bin/env bash
# get_subsample_data.sh

wget \
  -O "emp-paired-end-sequences/forward.fastq.gz" \
  "https://data.qiime2.org/2022.8/tutorials/atacama-soils/10p/forward.fastq.gz"

wget \
  -O "emp-paired-end-sequences/reverse.fastq.gz" \
  "https://data.qiime2.org/2022.8/tutorials/atacama-soils/10p/reverse.fastq.gz"

wget \
  -O "emp-paired-end-sequences/barcodes.fastq.gz" \
  "https://data.qiime2.org/2022.8/tutorials/atacama-soils/10p/barcodes.fastq.gz"
