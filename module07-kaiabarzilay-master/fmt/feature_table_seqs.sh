#!/usr/bin/env bash
# feature_table_seqs.sh

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
