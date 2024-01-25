#!/usr/bin/env bash
# clustalAlign.sh

clustalo --partition=express --job-name=clustalAlign --time=00:30:00 -N 1 -n 1 --output="batch-%x-%j.output" --infile={apoe_aa.fasta, -} --outfile={apoe_alignment.fasta, -}
