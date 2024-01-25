#!/usr/bin/env python3
# BioPython_genbank.py

from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

Entrez.email = "barzilay.k@northeastern.edu"

sequences_list = []

with Entrez.efetch(
	db="nucleotide", rettype="gb", retmode="text", id="515056"
) as handle:
	seq_record_1 = SeqIO.read(handle, "gb")
	
sequences_list.append(seq_record_1)

with Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text", id="J01673.1"
) as handle:
        seq_record_2 = SeqIO.read(handle, "gb")

sequences_list.append(seq_record_2)

for seq in sequences_list:
	print(seq.seq)

sequences_features = []

for seq_object in sequences_list:
	new_seq_feature = SeqFeature(location=FeatureLocation(seq_object.features.location),type=seq_object.features.type, strand=seq_object.features.strand)
	sequences_features.append(new_seq_feature)

for features in sequences_features:
	print("Type:" + features.type + ", Location:" + features.location + ", Strand:" + features.strand)








