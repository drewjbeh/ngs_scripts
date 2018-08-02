#! /usr/bin/env python3

import sys
import gzip

input_fq = sys.argv[1]
if "gz" in input_fq:
				input_fq = gzip.open(input_fq,'r')
else:
	input_fq = open(sys.argv[1],'r')

line_count = 0
read_count = 0
total_lengths = 0
min_length = 0
max_length = 0

for line in input_fq:
				line_count += 1
				if line_count == 2:
								min_length = len(line)
				if line_count % 4 == 2:
								read_count += 1
								total_lengths += len(line)
								if len(line) > max_length:
												max_length = len(line)
								if len(line) < min_length:
												min_length = len(line)

mean_length = total_lengths / read_count
print ("Mean read length: {}\n Max read length: {}\n Min read length: {}".format(mean_length, max_length, min_length))
