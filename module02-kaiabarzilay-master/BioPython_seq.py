#!/usr/bin/env python3
# BioPython_seq.py

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

# Create a SeqRecord object with seq, id, description and alphabet parameters

sequence_seq = Seq("AAAATGGGGGGGGGGGCCCCGTT", generic_dna)

from Bio.SeqRecord import SeqRecord

seq_1 = SeqRecord(sequence_seq, id = "#12345", description = "example 1")

print (seq_1.id)

# Write the object to a sequence file in GenBank format

SeqIO.write(seq_1, "BioPython_seq.gb", "genbank")
