#!/usr/bin/env python
# hamming.py

import sys

def hamming_distance(string1, string2):
	"""
	Caculate the hamming distance between two strings
	"""
	distance = 0
	L = len(string1)
	for i in range(L):
		if string1[i] != string2[i]:
			distance += 1
	return distance

if __name__ == "__main__":
	# Check to make sure there are at least two strings
	arg_count = len(sys.argv) - 1
	if arg_count < 2:
		raise Exception("This script requires 2 strings")
	string1 = sys.argv[1]
	string2 = sys.argv[2]
	distance = hamming_distance(string1, string2)
	print(string1 + "\t" + string2 + "\t" + str(distance))
