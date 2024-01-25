#!/usr/bin/env bash
# summarize.sh

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
