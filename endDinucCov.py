#!/usr/bin/env python3

# Calculate proportion of CCA vs CC at 3' ends of uniquely aligned tRNAseq reads using ratios of dinucleotides (i.e. CA vs CC)

import argparse
import pysam
from collections import defaultdict

parser = argparse.ArgumentParser(usage = "%(prog)s -b bamfile1, [bamfile2...]")

parser.add_argument('-b', metavar='bam files', required=True, dest='bams', help='Uniquely aligned reads in bam format.', nargs='*')
parser.add_argument('-o', metavar='output', required=True, dest='out', help='Output dinucleotide occurences.')
args = parser.parse_args()


with open(args.out, "w") as out:
	out.write("Alignment file\tDinucleotide\tProportion\n")
	for bam in args.bams:
		aln_count = 0
		bamname = bam.split("/")[-1].split(".")[0]
		dinuc_dict = defaultdict(int)
		currbam = pysam.AlignmentFile(bam, "rb")
		for read in currbam.fetch(until_eof=True):
			aln_count += 1
			dinuc = read.query_sequence[-2:]
			dinuc_dict[dinuc] += 1
		
		dinuc_dict_norm = {k: v / aln_count for k, v in dinuc_dict.items()}
		
		for k, v in dinuc_dict_norm.items():
			out.write(bamname + "\t" + k + "\t" + str(v) + "\n")