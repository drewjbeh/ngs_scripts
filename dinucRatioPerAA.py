#!/usr/bin/env python3

# Similar to endDinucCov.py but this calculates this ratio of CCA/CC (i.e. charged ratio, when > 1 more is charged, when < 1 more is uncharged)
# Calculates ratio for each cluster mapped to by mim-tRNAseq, and for each condition (e.g. WT vs mutant)
# Outputs table of results for plotting

import argparse
import pysam
from collections import defaultdict

parser = argparse.ArgumentParser(usage = "%(prog)s -s sampleFile -o outputFile")

parser.add_argument('-s', metavar='sample info', required=True, dest='samples', help='Sample information, tab-separated file with absolute bam location and group/condition')
parser.add_argument('-o', metavar='output', required=True, dest='out', help='Output file')
args = parser.parse_args()

with open(args.out, "w") as out:
	out.write("Condition\tCluster\tCharge ratio\n")
	with open(args.samples, "r") as samples:
		for sample in samples:
			sample = sample.strip()
			aln_count = 0
			bam = sample.split("\t")[0]
			bamname = bam.split("/")[-1].split(".")[0]
			condition = sample.split("\t")[1]
			charge_dict = defaultdict(lambda: defaultdict(int))
			currbam = pysam.AlignmentFile(bam, "rb")
			for read in currbam.fetch(until_eof=True):
				aln_count += 1
				dinuc = read.query_sequence[-2:]
				cluster = read.reference_name.split("-")[-4:]
				cluster = "-".join(cluster )
				if dinuc.upper() == 'CA':
					charge_dict[cluster]['CA'] += 1
				if dinuc.upper() == 'CC':
					charge_dict[cluster]['CC'] += 1
		
			charge_dict_norm = {
				outer_k: {
					inner_k: inner_v / aln_count 
					for inner_k, inner_v in outer_v.items()
				}
				for outer_k, outer_v in charge_dict.items()
			}

			charge_dict_norm = {
				outer_k: {
					'CR': outer_v['CA']/outer_v['CC']
				}
				for outer_k, outer_v in charge_dict_norm.items()
			} 

			for k, v in charge_dict_norm.items():
				out.write(condition + "\t" + k + "\t" + str(v['CR']) + "\n")
