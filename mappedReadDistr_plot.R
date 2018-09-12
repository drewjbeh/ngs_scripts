#!/usr/bin/env Rscript

# Plots mapped (aligned) read distributions from bam/sam files calculated using the script mappedReadDistr.sh 

library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
				stop("At least one argument must be supplied (input file).n", call.=FALSE)

} else if (length(args) > 0) {
				distr_plot_list = list()

				for (i in args ) {
								print(i)
								length_distr = read.table(i, header=FALSE)
								colnames(length_distr) = c("Frequency", "Length")

								distr_plot_list[[i]] = ggplot(length_distr, aes(x = Length, y = Frequency)) + 
								geom_bar(stat="identity") + scale_x_continuous(breaks = seq(10,80,5), limits = c(10,80)) + 
								ggtitle(strsplit(i,"_")[[1]][1])

				}

				length_distr_combinedPlot = do.call(grid.arrange, distr_plot_list)
				ggsave("mappedLengthDistribution.pdf",length_distr_combinedPlot, width = 12, height = 12, units = "in")
}
