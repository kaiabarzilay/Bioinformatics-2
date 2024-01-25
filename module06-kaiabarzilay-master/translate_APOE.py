#!/usr/bin/env python3
# translate_APOE.py
from Bio import SeqIO
from Bio.Seq import Seq

for seq_record in SeqIO.parse("APOE_refseq_transcript.fasta", "fasta"):
	amino_seq = seq_record.translate()
	SeqIO.write(amino_seq, "apoe_aa.fasta", "fasta")



