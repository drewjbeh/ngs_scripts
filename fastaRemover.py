#!/usr/bin/env python3

## Removes fasta records from a multi-fasta file based on IDs given in a second file

from Bio import SeqIO
import sys

inputFa = sys.argv[1]
inputIDs = open(sys.argv[2], 'r')

fasta = SeqIO.parse(inputFa,"fasta")
finalFasta = list()

IDs_list = list()

for line in inputIDs:
	line = line.rstrip('\n')
	IDs_list.append(line)

for seq in fasta:
	if seq.name not in IDs_list:
		finalFasta.append(seq)

SeqIO.write(finalFasta, inputFa.split(".")[0] + '_noDubious.fa', 'fasta')
