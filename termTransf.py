#!/usr/bin/env python3

import pysam
import sys, re
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

bam = sys.argv[1]
bam_name = bam.split("/")[-1].split(".")[0]
softclip_dict = defaultdict(int)
identity_dict = defaultdict(int)

bam_file = pysam.AlignmentFile(bam, "rb")
for read in bam_file.fetch(until_eof = True):
	cigar = read.cigarstring
	read_seq = read.query_sequence
	cigar_list = re.split('(.*?)([A-Za-z]|[\^][A-Za-z]+)', cigar)
	cigar_list = list(filter(None, cigar_list))
	if 'S' in cigar_list and cigar_list.index('S') == 1:
		softcliip_num = int(cigar_list[0])
		softclip_dict[softcliip_num] += 1
		identity_dict[read_seq[0]] += 1

softclip = plt.figure()
plt.bar(range(len(softclip_dict)), list(softclip_dict.values()), align = 'center')
plt.xticks(range(len(softclip_dict)), list(softclip_dict.keys()))
softclip.savefig(bam_name + "_term_distr.pdf", bbox_inches = "tight")

identity = plt.figure()
plt.bar(range(len(identity_dict)), list(identity_dict.values()))
plt.xticks(range(len(identity_dict)), list(identity_dict.keys()))
identity.savefig(bam_name + "_term_identity.pdf", bbox_inches = "tight")