#!/usr/bin/env python3

## Trim N's from 3' end of Ribosome profiling reads - useful for first RPy experiments (RPy_CHX_TIG_20180322) where reads had 
## Remaining 3' ends from PCR primer remanants

from Bio import SeqIO
import sys
import gzip
import subprocess

input = sys.argv[1]
output = sys.argv[2]

if ".gz" in input:
	subprocess.call(['gunzip', input])
	input = input.split(".gz")[0]

trimmed_reads = list()
n_count = 0

for record in SeqIO.parse(input, "fastq"):
	if not record.seq.startswith("NNNNNNNN"):
		while record.seq.endswith("N"):
			n_count += 1
			record = record[:-1]
		trimmed_reads.append(record)

count = SeqIO.write(trimmed_reads, output, "fastq")
print("Saved {} reads\n Removed {} N's".format(count, n_count))
