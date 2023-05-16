#!/usr/bin/env Rscript

#install.packages("ggrepel")
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))
library(ggrepel)

misincDEPlots <- function(pos, DEtable) {
  
  # create output folder
  sub_dir = paste0("extraPlots_",as.character(pos))
  output_dir = file.path(getwd(), sub_dir)
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    print("Dir already exists!")
  }
  
  # get control_cond and test_cond automatically from DESeq table name
  test_cond = unlist(strsplit(basename(DEtable), "vs"))[1]
  control_cond = unlist(strsplit(unlist(strsplit(basename(DEtable), "vs"))[2], "_diffexpr"))[1]
  pos = as.character(pos)
  mismatch_table = "mods/mismatchTable.csv"
  
  # read in mods context info and filter for pos
  modContext = read.table("mods/modContext.txt", header = TRUE)
  colnames(modContext)[1] = "isodecoder"
  modContext = modContext %>%
    select(-c(upstream, downstream)) %>%
    filter(canon_pos == pos) %>%
    select(-canon_pos)
  modContext$isodecoder = gsub(".*?_.*?_mito_tRNA-", "mito", modContext$isodecoder)
  modContext$isodecoder = gsub(".*?_.*?_tRNA-", "", modContext$isodecoder)
  modContext = modContext[!grepl("eColi", modContext$isodecoder),]
  
  # read in mods table and format/filter for pos
  mods = read.table(mismatch_table, header = TRUE)
  mods$isodecoder = gsub(".*?_.*?_mito_tRNA-", "mito", mods$isodecoder)
  mods$isodecoder = gsub(".*?_.*?_tRNA-", "", mods$isodecoder)
  mods = mods[!grepl("eColi", mods$isodecoder),]
  mods$proportion[is.na(mods$proportion)] = 0 # prevent losing NA values by setting to 0 here
  mods_pos = mods[mods$canon_pos == pos,]
  
  # aggregate by summing all 4 types and filter by conditions of interest
  mods_pos_agg = mods_pos %>% 
    select(-pos) %>% 
    group_by(isodecoder, condition, bam, cov, canon_pos) %>% 
    summarise(proportion = sum(proportion)) %>%
    filter(condition %in% c(control_cond, test_cond))
  
  # get mean over replicates
  mods_pos_agg = mods_pos_agg %>% group_by(isodecoder, condition, canon_pos) %>%
    summarise_each(funs(mean)) %>%
    select(-bam)
  
  # calculate misinc diff (= condition - control), mean cov, and create boolean for prop > 0.1 (thresh)
  mods_diff = mods_pos_agg %>% 
    group_by(isodecoder, canon_pos) %>% 
    filter(any(condition == control_cond)) %>% # keep isodecoders with data at control condition
    summarise(diff = proportion[condition == test_cond] - proportion[condition == control_cond], 
              cov = mean(cov),
              thresh = proportion[condition == control_cond] >= 0.1,
              control_misinc = proportion[condition == control_cond])
  
  # read abindance data from DESeq2 and format
  abundance = read.table(DEtable, header = TRUE, sep = ",")[,-c(1)]
  colnames(abundance)[1] = "isodecoder"
  abundance$isodecoder = gsub(".*?_.*?_mito_tRNA-", "mito", abundance$isodecoder)
  abundance$isodecoder = gsub(".*?_.*?_tRNA-", "", abundance$isodecoder)
  abundance = abundance[!grepl("eColi", abundance$isodecoder),]
  # new variables: sig is bool of significant DE, direction is direction of DE
  abundance$sig = ifelse(abundance$padj <= 0.05, "TRUE", "FALSE")
  abundance$direction = ifelse(abundance$padj <= 0.05, sign(abundance$log2FoldChange), 0)
  abundance$direction = as.character(abundance$direction)
  # merge with context info
  abundance = abundance %>%
    left_join(modContext, by = "isodecoder")
  
  # read in CCA data
  cca = read.table("CCAanalysis/CCAcounts.csv", header = T)
  colnames(cca)[colnames(cca) == "gene"] = "isodecoder"
  cca$isodecoder = gsub(".*?_.*?_mito_tRNA-", "mito", cca$isodecoder)
  cca$isodecoder = gsub(".*?_.*?_tRNA-", "", cca$isodecoder)
  cca = cca[!grepl("eColi", cca$isodecoder),]
  cca = cca %>%
    group_by(isodecoder, sample) %>%
    mutate(prop = count/sum(count))
  cca_avg = cca %>%
    group_by(isodecoder, end, condition) %>%
    summarise(prop_mean = mean(prop))
  cca$end = factor(cca$end, levels = c("CA", "CC", "C", "Absent"))
  cca_avg$end = factor(cca_avg$end, levels = c("CA", "CC", "C", "Absent"))
  
  cca_diff = cca_avg %>%
    group_by(isodecoder) %>%
    summarise(diff = prop_mean[end == "CA" & condition == test_cond] - prop_mean[end == "CA" & condition == control_cond],
      control_cca = prop_mean[end == "CA" & condition == control_cond]) %>%
    mutate(thresh = diff <= -0.1)
  
  # merge abundance with mods data
  merged_counts = abundance %>% select(-c(lfcSE, stat, pvalue, size)) %>%
    right_join(mods_pos_agg, by="isodecoder")
  merged_counts$condition = factor(merged_counts$condition, levels=c(`control_cond`, `test_cond`))
  
  # lookup table for classifying tRNAs as combinations of meeting misinc thesh (0,1) or being DE
  typeTable = data.frame(comb = c("TRUE_TRUE", "TRUE_FALSE", "FALSE_TRUE", "FALSE_FALSE"), 
                         type = c("both", "sig", "thresh", "none"))
  
  # merge abundance with mods diff data
  merged_diff  = abundance %>% select(-c(lfcSE, stat, pvalue, size)) %>%
    right_join(mods_diff, by = "isodecoder") %>%
    mutate(temp = paste(sig, thresh, sep = "_"))
  
  # match temp variable above to typeTable, remove temp
  merged_diff$type = typeTable$type[match(merged_diff$temp, typeTable$comb)]
  merged_diff = merged_diff %>% 
    select(-temp)
  merged_diff$type = factor(merged_diff$type, levels = c("none", "thresh", "sig", "both"))
  
  # merge abundance with cca diff data
  merged_cca = abundance %>% select(-c(lfcSE, stat, pvalue, size)) %>%
    left_join(cca_diff, by = "isodecoder") %>%
    mutate(temp = paste(sig, thresh, sep = "_"))
  
  merged_cca$type = typeTable$type[match(merged_cca$temp, typeTable$comb)]
  merged_cca = merged_cca %>% 
    select(-temp)
  merged_cca$type = factor(merged_cca$type, levels = c("none", "thresh", "sig", "both"))
  
  # set width of plot based on facet number (unique values in identity column)
  facets = unique(merged_diff$identity)
  facets = length(facets[!is.na(facets)])
  width = ifelse(facets == 1, 1, 2)
  height = ceiling(facets / 2)
  # plot misinc diff over abundance fold change with context info
  diffFoldChange = ggplot(subset(merged_diff, !is.na(padj)), aes(x = diff, y = log2FoldChange)) + 
    geom_point(aes(color = type, size = control_misinc)) + #size = sig)) + 
    geom_text_repel(data = subset(merged_diff, !is.na(padj) & sig == TRUE ), 
              aes(label = isodecoder),
              size = 2.5, 
              hjust = "inward",
              vjust = "inward") +
    theme_bw() +
    facet_wrap(~identity, ncol = 2) +
    ylab(paste("Mature tRNA log2FC (", test_cond, "/", control_cond, ")", sep = )) +
    xlab(paste("Misincorporation difference (", test_cond, " - ", control_cond, ")", sep = )) +
    scale_color_manual(labels = c(none = "None", thresh = "Control misinc. >= 0.1", sig = "Diff. expressed (p-adj <= 0.05)", both = "Both"), 
                       values = list(none = "#bababa", thresh = "#404040", sig = "#f4a582", both = "#ca0020")) +
    scale_size_binned(range = c(0.05, 3), name = "Control misinc. proportion") +
    #scale_size_manual(labels = c(none = "None", thresh = "Control misinc. >= 0.1", sig = "Diff. expressed (p-adj <= 0.05)", both = "Both"), values = c(1, 2, 2, 2), guide = "none") +
    labs(color = "Type", size = "Misinc. proportion\nin control") +
    ggtitle(paste0("position ", pos)) + 
    xlim(c(-1, max(merged_diff$diff)))
  ggsave(paste0(output_dir, "/", control_cond, 'vs', test_cond, "_pos", pos, "_misincDiffLFC.pdf"), diffFoldChange, width = (7 * width) + 1, height = 6 * height, device = "pdf")
  
  # plot absolute misinc props with change in abundance
  pd <- position_dodge(0.4)
  misincSign = ggplot(subset(merged_counts, !is.na(padj)), aes(x = condition, y = proportion, group = isodecoder)) + 
    geom_line(data = subset(merged_counts, padj <= 0.05), position = pd, linetype = "dashed", alpha = 0.4, color = "#C14633") + 
    geom_point(position = pd, aes(color = direction, shape = direction, size = direction)) +
    theme_bw() +
    scale_color_manual('Differential expression\n(adj-p <= 0.05)', labels = c("Down", "None", "Up"), values = c("#f28f3b", "darkgrey", "#4daf4a")) + 
    scale_shape_manual('Differential expression\n(adj-p <= 0.05)', labels = c("Down", "None", "Up"), values = c(17, 19, 17)) +
    scale_size_manual('Differential expression\n(adj-p <= 0.05)', labels = c("Down", "None", "Up"), values = c(3, 1, 3)) +
    geom_text_repel(data = subset(merged_counts, direction !=0 & condition == test_cond), 
              aes(label = isodecoder),
              size = 3.5,
              hjust = -0.3,
              position = pd) +
    facet_wrap(~identity, ncol = 2) +
    ylab(paste0("Misincorporation proportion at position ", pos)) +
    theme(axis.title.x = element_blank()) +
    ggtitle(paste0("position ", pos)) +
    ylim(c(0,1))
  ggsave(paste0(output_dir, "/", control_cond, 'vs', test_cond, "_pos", pos, "_misincPropsDE.pdf", sep = ""), misincSign, width = 14, height = 10, device = "pdf")
  
  # plot cca end type faceted by condition
  ccaEndType = ggplot(cca_avg, aes(x = end, y = prop_mean, fill = end)) +
    #geom_jitter(data = cca, aes(x=end, y=prop, color = sample), show.legend = FALSE, alpha = 0.6, width = 0.3) +
    geom_jitter(data = cca, aes(x=end, y=prop), color = "#bdbdbd", show.legend = FALSE, alpha = 0.6, width = 0.3) +
    geom_violin(alpha = 0.6, show.legend = FALSE, outlier.shape = NA, scale = "width", draw_quantiles = c(0.5)) + 
    geom_text_repel(data = subset(cca_avg, (end == 'CA' & prop_mean < 0.80) | (end == 'CC' & prop_mean > 0.2)), 
              aes(label = isodecoder),
              size = 3.5,
              hjust = "inward",
              vjust = "inward") +
    facet_grid(~condition) +
    #scale_color_manual(values = c("#bdbdbd", "#636363","#bdbdbd", "#636363")) +
    scale_fill_manual(name = "", values = alpha(c(CA = "#F0F9ED", CC = "#427AA1", C = "#0D4282", Absent = "#133C55"), 0.8)) +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    ylab("Mean end type proportion")
  ggsave(paste0(output_dir, "/", control_cond, 'vs', test_cond, "_pos", pos, "_CCAendProps.pdf", sep = ""), ccaEndType, width = 10, height = 6, device = "pdf")
  
  # plot cca (only) end proportion, x = condition
  cca_avg_filt = cca_avg[(cca_avg$end == "CA") & !grepl("mito", cca_avg$isodecoder),] %>%
                    merge(cca_diff, by = "isodecoder") %>%
                    mutate(label = ifelse(thresh == TRUE,
                                          isodecoder, ""))
  cca_filt = cca[(cca$end == "CA") & !grepl("mito", cca$isodecoder),] 

  ccaEndByCond = ggplot(cca_filt, aes(x = condition, y = prop, fill = condition)) +
    #geom_jitter(data = cca, aes(x=end, y=prop, color = sample), show.legend = FALSE, alpha = 0.6, width = 0.3) +
    geom_jitter(data = cca_filt, aes(x=condition, y=prop, color = sample), show.legend = FALSE, alpha = 0.6, width = 0.1) +
    geom_violin(alpha = 0.6, show.legend = FALSE, outlier.shape = NA, scale = "width", draw_quantiles = c(0.5)) + 
    geom_text_repel(data = cca_avg_filt,
                  aes(label = label, x = condition, y = prop_mean), 
                  max.overlaps = 100,
                  #box.padding = 1,
                  #show.legend = FALSE,
                  size = 2,
                  hjust = "inward",
                  vjust = "inward") +
    scale_color_manual(values = c("#bdbdbd", "#bdbdbd", "#636363", "#636363")) +
    scale_fill_manual(name = "", values = alpha(c("#C84630", "#427AA1"), 0.8)) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.title.x = element_blank()) +
    ylab("Mean 3'-CCA proportion")
  ggsave(paste0(output_dir, "/", control_cond, 'vs', test_cond, "_pos", pos, "_CCAbyCondition.pdf", sep = ""), ccaEndByCond, width = 6, height = 4, device = "pdf")
  

  # plot cca diff with abundance change info
  ccaDiffFoldChange = ggplot(subset(merged_cca, !is.na(padj)), aes(x = diff, y = log2FoldChange)) + 
    geom_point(aes(color = type, size = control_cca)) + 
    geom_text_repel(data = subset(merged_cca, !is.na(padj) & thresh == TRUE ), 
              aes(label = isodecoder),
              size = 3.5, 
              hjust = "inward",
              vjust = "inward") +
    theme_bw() +
    #facet_wrap(~identity, ncol = 2) +
    ylab(paste("Mature tRNA log2FC (", test_cond, "/", control_cond, ")", sep = )) +
    xlab(paste("Difference in 3'-CCA proportion (", test_cond, " - ", control_cond, ")", sep = )) +
    scale_color_manual(labels = c(none = "None", thresh = "3'-CCA diff >= 0.1", sig = "Diff. expressed (p-adj <= 0.05)", both = "Both"), 
                       values = list(none = "#bababa", thresh = "#f4a582", sig = "#404040", both = "#ca0020")) +
    scale_size_binned(range = c(0.05, 3), breaks = c(0,0.75,1), name = "Control 3'-CCA proportion") +
    #scale_size_manual("3-CCA difference >= 0.1", labels = c("FALSE" = "False", "TRUE" = "True"), values = c(1, 2), guide = "none") +
    ggtitle(paste0("position ", pos)) +
    xlim(c(-1, max(merged_cca$diff)))
  ggsave(paste0(output_dir, "/", control_cond, 'vs', test_cond, "_pos", pos, "_CCADiffLFC.pdf", sep = ""), ccaDiffFoldChange, width = 10, heigh = 6, device = "pdf")
  
  #readout tables
  write.csv(merged_diff,paste(output_dir,"merged_mods_diff.csv", sep = "/"))
  write.csv(merged_counts,paste(output_dir,"merged_mods_counts.csv", sep = "/"))
  write.csv(cca_avg,paste(output_dir,"cca_avg.csv", sep = "/"))
  write.csv(merged_cca,paste(output_dir,"merged_cca.csv", sep = "/"))
  
}

# args!
parser <- ArgumentParser(description = "Plot changes in misincorporation and transcript abundance at specific tRNA position between two conditions")
parser$add_argument("pos", nargs=1, help="tRNA position to analyse")
parser$add_argument("DEseq", nargs=1, help="Path to mim-tRNAseq isodecoder DESeq results for the comparison. Can be cytosolic/mitochondrial results file.")
args <- parser$parse_args()
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
pos = as.integer(args$pos)
DEtable = args$DEseq

misincDEPlots(pos, DEtable)
