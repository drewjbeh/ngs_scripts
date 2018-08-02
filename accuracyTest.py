#! /usr/bin/env python3

# Calculations of accuracy (true positive, false positive) from simulated data mapped with tRNAseq pipeline

import argparse
import pysam

parser = argparse.ArgumentParser(description = 'Test accuracy of mapped simulated data by outputting important false and true negative and positives, and plots.')
parser.add_argument('-b','--bam', metavar = 'input bam files', required = True, dest = 'bams', nargs = '*')
parser.add_argument('-c','--clusterInfo', metavar='cluster information', required = False, dest= 'clusters', \
	help = 'Cluster information file, tab-delimited with two columns - gene name and cluster number.')
args = parser.parse_args()

def calcScores():

	if args.clusters:
		clusters = dict()
		with open(args.clusters,"r") as clusterInfo:
			for line in clusterInfo:
				line = line.strip()
				clusters[line.split("\t")[0]] = line.split("\t")[1]

	for bam in args.bams:
		total = 0
		tp = 0
		fp = 0
		print(bam)
		bam_file = pysam.AlignmentFile(bam, "rb")
		for read in bam_file.fetch(until_eof=True):
			total += 1
			query = read.query_name.split("/")[1]
			reference = read.reference_name
			if not query == reference:
				if args.clusters and clusters[query] == clusters[reference]:
					tp += 1
				else:
					fp += 1
			elif query == reference:
				tp += 1
		print("True positives: {} ({:.0%})".format(tp, tp/total))
		print("False positives: {} ({:.0%})".format(fp, fp/total))
		print("Total: {}".format(total))

calcScores()