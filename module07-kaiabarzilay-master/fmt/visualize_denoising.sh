#!/usr/bin/env bash
# visualize_denoising.sh

qiime metadata tabulate \
  --m-input-file stats-1.qza \
  --o-visualization denoising-stats-1.qzv
qiime metadata tabulate \
  --m-input-file stats-2.qza \
  --o-visualization denoising-stats-2.qzv
