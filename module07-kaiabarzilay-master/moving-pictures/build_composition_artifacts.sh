#!/usr/bin/env bash
# build_composition_artifacts.sh

qiime composition add-pseudocount \
  --i-table gut-table.qza \
  --o-composition-table comp-gut-table.qza
