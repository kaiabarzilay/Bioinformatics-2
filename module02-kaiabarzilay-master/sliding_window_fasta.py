#!/usr/bin/env python3
# sliding.py

from Bio import SeqIO
import sys


def get_kmers(k, string):
        """
        Returns a list of all kmers in a given string
        """
        kmers = []
        end = len(string) - k + 1
        for start in range(0, end):
                kmers.append(string[start:start + k])
        return kmers

def gc_content(dna):
	''' returns the fraction of GC in a string '''
	dna = dna.lower()
	gc = 0
	for nucleotide in dna:
		if nucleotide in ['g', 'c']:
			gc += 1
	return gc/len(dna)

if __name__ == "__main__":
        # Check to make sure there are at least two arguments
	arg_count = len(sys.argv) - 1
	if arg_count < 2:
		raise Exception("This script requires 2 arguments: a kmer size and then a fasta file")
	k = sys.argv[1]
	seq = sys.argv[2]
	kmers = get_kmers(int(k),seq)

	# Use SeqIO to open fasta file
	for record in SeqIO.parse(sys.argv[2], "fasta"):
		print(">{}".format(record.description))
		for kmer in get_kmers(int(sys.argv[1]), record.seq):
			print("{}\t{:.2f}".format(kmer, gc_content(kmer)))

