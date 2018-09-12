#!/usr/bin/env bash

# Calculate frequency table of mapped read length
# useful to plot in R using mappedReadDistr_plot.R
# Usage: mappedReadDistr.sh file1 [file2 file3 ...]

for FILE in "$@"
do
				FN_without_extension=${*/%FILE%.*}
				FN_without_extension=${FILE##*/}
				samtools view $FILE | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c > $FN_without_extension.distr
done	
