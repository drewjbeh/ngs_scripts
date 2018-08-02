#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
library(scales)

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)

} else if (length(args) > 0) {
  length_plot_list = list()
  frame_plot_list = list()
  
  for (i in args) {
    length_frame_distr = read.csv(i, header = TRUE)
    colnames(length_frame_distr) = c('Length',"Mismatch", "Frame",'Count')
    length_frame_distr$Frame = factor(length_frame_distr$Frame, levels = c("2","1","0"))
    
    length_plot_list[[i]] = ggplot(length_frame_distr, aes(x = Length, y = Count, group = Mismatch, fill = Mismatch)) + 
      geom_bar(stat = "identity") + scale_fill_discrete(drop=FALSE, name="5' mismatch") + 
      scale_x_continuous(breaks = seq(10,40,2), limits = c(10,40)) + theme_gray() + ggtitle(strsplit(i,"_")[[1]][1])
    
    frame_plot_list[[i]] = ggplot(length_frame_distr, aes(x = Length, y = Count, group = Frame, fill = Frame)) +
      facet_grid(~ Mismatch, scales = "free") +
      geom_bar(stat = "identity") + scale_fill_discrete(drop=FALSE) + 
      scale_x_continuous(breaks = seq(10,40,3), limits = c(10,40)) + theme_gray() + ggtitle(strsplit(i,"_")[[1]][1]) + 
      scale_fill_manual(values = c('2' = '#619CFF', '1' = '#F8766D', '0' = '#00BA38'))
  } 
  
  length_distr_combinedPlot = do.call(grid.arrange, length_plot_list)
  frame_distr_combinedPlot = do.call(grid.arrange, frame_plot_list)
  ggsave("LengthDistribution.pdf",length_distr_combinedPlot, width = 12, height = 12, units = "in")
  ggsave("FrameDistribution.pdf", frame_distr_combinedPlot, width = 12, height = 12, units = "in")
}
