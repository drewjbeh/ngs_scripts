#!/usr/bin/env python3

## Reads bam files (typically for ribo-seq) and analyses various attributes such as length, 5' mismatch nucleotide and frame
## Calls LengthandFramePlot.R -- make sure it is in the same directory as this script when running

import sys
import numpy as np
import pandas as pd
from collections import defaultdict
import pysam
import argparse
import re
import subprocess
import gffutils

parser = argparse.ArgumentParser(description = 'Analyse Ribo-seq bam files and extract mapped read characteristics')
parser.add_argument('-b','--bam', metavar = 'input bam files', required = True, dest = 'bams', nargs = '*')
parser.add_argument('-n', '--name', metavar = 'name for file prefix', required = True, dest = 'name')
parser.add_argument('--read-min', metavar = 'min read length', required = False, dest = 'min', help = 'Min read length to include in plots')
parser.add_argument('--read-max', metavar = 'max read length', required = False, dest = 'max', help = 'Max read length to include in plots')
parser.add_argument('--mismatch', dest = 'mismatch', required = False, action = 'store_true',\
	help = "Plot 5' mismatch data")
#parser.set_defaults(min = 10, max = 45)
args = parser.parse_args()

mapped_count = 0
results_list = list()

for bam in args.bams:
	seq_dict = defaultdict(list)
	print(bam)
	bam_file = pysam.AlignmentFile(bam,"rb")
	for read in bam_file.fetch(until_eof=True):
		if not read.is_unmapped:
			mapped_count += 1
			seq_dict['length'].append(read.query_length)
			seq_dict["5'_mis"].append(bool(re.search('^0[A:Z].',read.get_tag("MD"))))
			seq_dict['frame'].append(read.reference_start % 3)

	seq_df = pd.DataFrame(seq_dict)
	
	frame_size_distr = seq_df.groupby(['length',"5'_mis","frame"])['length'].count()
	frame_size_distr = pd.DataFrame(frame_size_distr)
	frame_size_distr.columns = ['count']

	filename = bam.split("/")[-1].split(".")[0]
	frame_size_distr.to_csv(filename + "_LengthFrame.csv")
	results_list.append(filename + "_LengthFrame.csv")

# set default min and max read lengths based on those in data
min_length = min(frame_size_distr.index.get_level_values(0))
max_length = max(frame_size_distr.index.get_level_values(0))
parser.set_defaults(min = min_length, max = max_length)
args = parser.parse_args()

command = "Rscript"
path = "./LengthandFramePlot.R"
cmd = [command, path, str(args.min), str(args.max), str(args.mismatch), str(args.name)] + results_list
x = subprocess.check_output(cmd, universal_newlines=True)
	#size_distr.to_csv(filename + "_LengthDistr.csv")
	#frame_distr.to_csv(filename + "_FrameDistr.csv")



