#!/usr/bin/env python3

import sys

input = open(sys.argv[1],'r')
output = open(sys.argv[2], 'w')

for line in input:
				line = line.strip()
				if line.split('\t')[2] == 'gene':
								id = line.split('\t')[8].split(';')[0].split('"')[1]
								name = line.split('\t')[8].split(';')[1].split('"')[1]
								out = line.split('\t')[0] + '\t' + line.split('\t')[1] + '\t' + line.split('\t')[2] + '\t' + line.split('\t')[3] + '\t' + line.split('\t')[4] + '\t' + '1000' + '\t' + line.split('\t')[6] + '\t' + line.split('\t')[7] + '\t' + 'ID=Gene:' + id + ';' + "Name=" + name + '\n'
								output.write(out)
				else:
								id = line.split('\t')[8].split(';')[1].split('"')[1]
								parent = line.split('\t')[8].split(';')[0].split('"')[1]
								out = line.split('\t')[0] + '\t' + line.split('\t')[1] + '\t' + line.split('\t')[2] + '\t' + line.split('\t')[3] + '\t' + line.split('\t')[4] + '\t' + '1000' + '\t' + line.split('\t')[6] + '\t' + line.split('\t')[7]     + '\t' + 'ID=' + id + ';' + "Parent=" + parent + '\n'
								output.write(out)
input.close()
output.close()
