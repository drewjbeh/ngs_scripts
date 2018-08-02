#!/usr/bin/env python3

## Extract records from gff matching those in match list - useful for aligning and filtering sequening reads to various noncoding loci

import sys

input = open(sys.argv[1],'r')
output = open(sys.argv[2],'w')

match = ["ncRNA","pseudogene","rRNA","snoRNA","snRNA","transposable_element_gene","tRNA"]
for line in input:
				line = line.strip()
				if not line.startswith("#"):
								if line.split("\t")[2] in match:
												output.write(line + "\n")

input.close()
output.close()
