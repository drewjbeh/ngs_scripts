#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
library(scales)

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)", call=FALSE)

} else if (length(args) > 0) {
  length_plot_list = list()
  frame_plot_list = list()

  min_length = args[1]
  max_length = args[2]
  mismatch = args[3]
  name = args[4]
  sample_num = length(5:length(args))

  length_frame_distr = data.frame(Length=integer(),Mismatch=logical(), Frame=integer(),Count=integer(), bam=character())
  
  for (i in 5:length(args)) {
    print(args[i])
    file = args[i]
    temp = read.csv(file, header = TRUE)
    file_name = strsplit(file,"transcriptome")[[1]][1]
    temp$bam = file_name
    colnames(temp) = c("Length","Mismatch", "Frame","Count", "bam")
    length_frame_distr = rbind(length_frame_distr, temp)
  }

  length_frame_distr$Frame = factor(length_frame_distr$Frame, levels = c("2","1","0"))

  length_plot = ggplot(length_frame_distr, aes(x = Length, y = Count, group = Mismatch, fill = Mismatch)) + 
    geom_bar(stat = "identity") +
    facet_wrap(~ bam, scales = "free", nrow = ceiling(sample_num/3)) +
    scale_fill_discrete(drop=FALSE, name="5' mismatch") + 
    scale_x_continuous(breaks = seq(as.numeric(min_length), as.numeric(max_length),3), limits = c(as.numeric(min_length), as.numeric(max_length))) + 
    theme_bw()
    
  if (mismatch == 'True'){
    frame_plot = ggplot(length_frame_distr, aes(x = Length, y = Count, group = Frame, fill = Frame)) +
    facet_grid(bam~Mismatch, scales = "free") +
    geom_bar(stat = "identity") + 
    scale_fill_discrete(drop=FALSE) + 
    scale_x_continuous(breaks = seq(as.numeric(min_length), as.numeric(max_length),3), limits = c(as.numeric(min_length), as.numeric(max_length))) + 
    theme_bw() + 
    scale_fill_manual(values = c('2' = '#619CFF', '1' = '#F8766D', '0' = '#00BA38'))
  } else if (mismatch == 'False'){
    frame_plot = ggplot(length_frame_distr, aes(x = Length, y = Count, group = Frame, fill = Frame)) +
    facet_wrap(~bam, scales = "free", nrow = ceiling(sample_num/3)) +
    geom_bar(stat = "identity") + 
    scale_fill_discrete(drop=FALSE) + 
    scale_x_continuous(breaks = seq(as.numeric(min_length), as.numeric(max_length),3), limits = c(as.numeric(min_length), as.numeric(max_length))) + 
    theme_bw() + 
    scale_fill_manual(values = c('2' = '#619CFF', '1' = '#F8766D', '0' = '#00BA38'))
  }

  ggsave(paste(name, "LengthDistribution.pdf", sep="_"),length_plot, width = 12, height = 16, units = "in")
  ggsave(paste(name, "FrameDistribution.pdf", sep="_"), frame_plot, width = 12, height = 16, units = "in")
}
