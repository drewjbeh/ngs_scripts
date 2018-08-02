#! /usr/bin/env python3

# Calculate occurence of dinucleotides at end of trimmed reads

import argparse
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser(usage = "%(prog)s -f fastafile1, [fastafile2...]")

parser.add_argument('-f', metavar='trimmed fastq files', required=True, dest='fq', help='Trimmed fastq files to analyze.', nargs='*')
args = parser.parse_args()

with open("dinuc_occurence.txt","w") as out:
	out.write("lib\tdinucleotide\tproportion\n")
	for fq in args.fq:
		linecount = 0
		readcount = 0
		Ncount = 0
		dinucl_dict = {'AA':0, 'AC':0, 'AG':0, 'AT':0, 'CA':0, 'CC':0, 'CT':0, 'CG':0, 'GA':0, 'GC':0, 'GG':0, 'GT':0, 'TA':0, 'TC':0, 'TG':0, 'TT':0}
		with gzip.open(fq, 'rt') as file:
			for line in file:
				line = line.strip()
				linecount += 1
				if linecount % 4 == 2:
					dinucl = line[-2:]
					if not "N" in dinucl:
						readcount += 1
						dinucl_dict[dinucl] += 1
					elif "N" in dinucl:
						Ncount += 1

		print("dinucleotides containing 'N' and not analysed: " + str(Ncount))

		dinucl_dict_norm = {k: v / readcount for k, v in dinucl_dict.items()}
		lib = fq.split("/")[-1].split(".")[0]

		for k, v in dinucl_dict_norm.items():
			out.write(lib + "\t" + k + "\t" + str(v) + "\n")


