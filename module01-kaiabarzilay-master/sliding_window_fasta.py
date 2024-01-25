#!/usr/bin/env python
# sliding.py

import re
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

def gc_content(kmers):
	gc_count = []
	for start in range(0, len(kmers)):
		string = kmers[start].lower()
		gc = 0
		for nucleotide in string:
			if nucleotide in ['g', 'c']:
		# Added due to python3 not working on my computer momentarly, 
                # because in python2 integers don't automatically convert to floats. 
                # Means any decimal will be rounded to zero		
				gc += 1
		gc = float(gc)
		string_length = float(len(string))
		gc_count.append(round(gc/string_length, 2))
	return gc_count

if __name__ == "__main__":
	# Check to make sure there are at least two arguments
	arg_count = len(sys.argv) - 1
	if arg_count < 2:
		raise Exception("This script requires 2 arguments: a kmer size and then a fasta file")
	elif int(sys.argv[1]) <= 0:
        	sys.exit("Value of K must be greater then Zero")
	k = int(sys.argv[1])
	file_fasta = sys.argv[2]
	kmer = ""

with open(file_fasta) as genome:
	"""
	Opens file and gets the sequence in a string, printing the header.
	"""
        for line in genome:
            line = line.strip()
            if len(line) < 1:
                continue
            elif re.match('^>', line):
                print(line)
                kmer = ""
            else:
                kmer += line
	kmers = get_kmers(k, kmer)
        content = gc_content(kmers)
        for i in range(len(kmers)):
		print("{}\t{:.2f}".format(kmers[i], content[i]))


