#!/usr/bin/env bash
# export_results.sh

qiime tools export \
  --input-path feature-table.qza \
  --output-path exported-feature-table
