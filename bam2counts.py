#!/usr/bin/env python3

# Counts reads mapping to features in input bam files

import sys
import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description = 'Count mapped reads per feature in bam files, producing count table')
parser.add_argument('-b','--bam', metavar = 'input bam files', required = True, dest = 'bams', nargs = '*')
parser.add_argument('-o', '--out', metavar = 'output file prefix', required = False, dest = 'out')
parser.set_defaults(out = '')
args = parser.parse_args()

if args.out != '':
	output = args.out + "_countsPerFeature.csv"
else:
	output = "countsPerFeature.csv"

with open(output, "w") as out:
	out.write("Library,Target,Count\n")
	for bam in args.bams:
		mapped_count = 0
		unmapped_count = 0
		total = 0
		count_dict = defaultdict(int)
		library = bam.split("/")[-1].split(".")[0]
		bam_file = pysam.AlignmentFile(bam,"rb")
		for read in bam_file.fetch(until_eof=True):
			if not read.is_unmapped:
				mapped_count += 1
				total += 1
				count_dict[read.reference_name] += 1
			else:
				total += 1
				unmapped_count += 1
		for feature in count_dict:
			out.write(library + "," + feature + "," + str(count_dict[feature]) + "\n")
		print('{} processed: {} total reads, {} reads mapped to {} features'.format(bam, str(total), str(mapped_count), len(count_dict)))

