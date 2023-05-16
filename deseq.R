library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(calibrate)
library(plyr)
library(grid)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

lastlevel = function(f, control) {
  if (!is.factor(f)) stop("input for contrast combinations not a factor")
  orig_levels = levels(f)
  if (! control %in% orig_levels) stop("Control must be a level of data - potentially removed due to single replicate. Control condition must have replicates!")
  new_levels = c(setdiff(orig_levels, control), control)
  factor(f, levels = new_levels)
}

control_condition = "kips_nt_t5"
p_adj = 0.05
outdir = "./"
counts = read.csv("../4-counts/ProteinCoding.counts.txt", header=TRUE, sep="\t", comment.char = "#")
counts = counts[,-c(2,3,4,5,6)]
colnames(counts) = c("Geneid","kips_nt_t5_1", "kips_nt_t5_2","kips_gigyf2_t5_1","kips_gigyf2_t5_2")
rownames(counts) = counts$Geneid
counts = counts[,-1]

#coldata = data.frame(row.names = colnames(counts), condition = rep(c("dox", "nodox"), 2))
coldata = data.frame(row.names = colnames(counts), condition = c("kips_nt_t5", "kips_nt_t5","kips_gigyf2_t5","kips_gigyf2_t5"))
                          
dds = DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds = DESeq(dds)
baseMeanPerLvl = as.data.frame(sapply(levels(dds$condition), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl] )))
count_df = as.data.frame(counts(dds, normalized=TRUE))

# plots dispersion, save to out
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()
  
# shrink lfc estimates...
resLFC <- lfcShrink(dds, coef = 2, type="apeglm")
png("qc-MAplot.png", 1000, 1000, pointsize = 20)
plotMA(resLFC, ylim=c(-2,2))
dev.off()

#vsd transform
vsd = varianceStabilizingTransformation(dds, blind=FALSE)
assay_vsd = as.data.frame(assay(vsd)) %>% rownames_to_column(var = "Gene")
write.csv(assay_vsd, file = paste("vst-transformedCounts.csv", sep = "/"), row.names = FALSE)

# plot PCA and save to out
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
ggsave("qc-pca.pdf", height = 7, width = 8)
  
# create ordered levels so that the DE comparison is always condition vs control
# needs to be manually edited for each experiment with different cell types or DE analyses
ordered_levels = levels(lastlevel(unique(dds$condition), control_condition))
combinations = combn(ordered_levels, 2, simplify=FALSE)
DEcounts_list = list()

for (i in 1:length(combinations)) {
  res = results(dds, contrast=c("condition",as.vector(combinations[[i]])))
  res = res[order(res$padj), ]
  res$rn = rownames(res)
  normcount_df = as.data.frame(counts(dds, normalized=TRUE))
  normcount_df$rn = rownames(normcount_df)
  resdata = join_all(list(as.data.frame(res), normcount_df), by="rn", type = 'left')
  names(resdata)[7] = "Gene"
  col_idx = grep("Gene", names(resdata))
  resdata = resdata[, c(col_idx, (1:ncol(resdata))[-col_idx])]
    
  # Count plots
  # add significance to baseMean matrices for current contrast
  baseMeanPerLvl$sig = res[rownames(baseMeanPerLvl),'padj'] < p_adj
  baseMeanPerLvl$sig = !is.na(baseMeanPerLvl$sig) & baseMeanPerLvl$sig # deal with NAs turning them into FALSE
    
  # add direction of DE to baseMean
  baseMeanPerLvl$direction = sign(res[rownames(baseMeanPerLvl), 'log2FoldChange'])
  baseMeanPerLvl[is.na(baseMeanPerLvl$direction),'direction'] = 0
  
  baseMeanPerLvl$comb = paste(baseMeanPerLvl$sig, baseMeanPerLvl$direction, sep='')
  baseMeanPerLvl$comb[which(baseMeanPerLvl$comb %in% c('FALSE1','FALSE0','FALSE-1'))] = "FALSE"
  
  cor = format(cor(baseMeanPerLvl[,combinations[[i]][1]], baseMeanPerLvl[,combinations[[i]][2]]), digits = 2)
  
  # Convert x and y variables into symbols  to handle str and numeric variable names
  x = sym(combinations[[i]][2])
  y = sym(combinations[[i]][1])
  count_data = subset(baseMeanPerLvl, select = c(combinations[[i]], "sig", "comb"))
  
  diff_dot = ggplot(count_data, aes_string(x, y)) +
    geom_point(aes(color = comb, shape = comb, size = sig)) +
    scale_x_log10() + scale_y_log10() + 
    #geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, color = '#3182bd', fill = 'grey') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = '#3182bd', alpha = 0.8) + 
    scale_color_manual('Differential expression', labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
    scale_shape_manual('Differential expression', labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
    scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() +
    labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
    annotate("label", 0, Inf, hjust = 0, vjust = 1, label = paste("italic(r) == ", cor), parse = TRUE)
    
  ggsave(paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), diff_dot, height = 5, width = 8)
    
  write.csv(resdata, file=paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"))
  
  if (control_condition %in% combinations[[i]]) {
    temp_lfc = as.data.frame(res[,c("log2FoldChange","padj")])
    temp_lfc = subset(temp_lfc, padj <= p_adj)
    colnames(temp_lfc) = c(paste(paste(combinations[[i]], collapse = "vs"), "_l2FC", sep = ""), paste(paste(combinations[[i]], collapse = "vs"), "_padj", sep = ""))
    temp_lfc = tibble::rownames_to_column(temp_lfc, var = "gene")
    temp_lfc = temp_lfc[!grepl("Escherichia_coli", temp_lfc$gene),]
    DEcounts_list[[paste(combinations[[i]], collapse = "vs")]] = temp_lfc
  }
}

# Build combined filtered tables for DE heatmaps and scale counts
comb = DEcounts_list %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene"), .)
basemean = as.data.frame(res) %>% dplyr::select(baseMean) %>% tibble::rownames_to_column(var = "gene")
normcounts = count_df %>% tibble::rownames_to_column(var = "gene")
comb = list(comb, basemean, normcounts) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="gene"), .)

# scaled Z-score expression matrix
scaled_counts = as.matrix(t(scale(t(comb[,!grepl("gene|l2FC|padj|baseMean", colnames(comb))]))), rownames.force = TRUE)
scaled_counts[is.na(scaled_counts)] = 0
rownames(scaled_counts) = comb$gene

# Heatmaps
# note that all tables and plots only built if DE tables are not empty
col_fun = colorRamp2(c(-3, 0, 3), c("#762a83", "#f7f7f7", "#1b7837"))

# Annotation for baseMean counts
if (nrow(comb) != 0) {
  baseMeanAnno_iso = rowAnnotation(base_mean = anno_lines(comb$baseMean, 
                                                          gp = gpar(lwd = 1.5, col = "#084081")), 
                                   annotation_name_rot = 90,
                                   width = unit("0.8", "cm"),
                                   height = unit(20, "cm"))
  
  hm_gene = Heatmap(scaled_counts,
                   col = col_fun,
                   row_title = paste('n = ', nrow(scaled_counts), sep = ""),
                   show_row_names = FALSE,
                   width = unit((length(ordered_levels)*0.7) + 2, "cm"),
                   height = unit(20, "cm"),
                   border = "gray20",
                   cluster_columns = TRUE,
                   heatmap_legend_param = list(
                     title = paste0("Scaled expression (Z score)\nlog2 fold-change (adj-p <= ", as.character(p_adj), ")")
                   )
  )
}

for (i in seq_len(length(combinations))){
  if (control_condition %in% combinations[[i]]) {
    lfc = paste(paste(combinations[[i]], collapse="vs"), "l2FC", sep="_")
    if (nrow(comb) != 0) {
      lfc_gene = Heatmap(comb[lfc],
                        col = col_fun,
                        width = unit(0.5, "cm"),
                        height = unit(20, "cm"),
                        na_col = "white",
                        border = "gray20",
                        show_heatmap_legend = FALSE)
      hm_gene = hm_gene + lfc_gene
    }
  }
}

# calculate plot width based on sample numbers
plot_width = (length(ordered_levels)*0.8 + 4) + (length(ordered_levels)-1)*1.3 + 0.8
plot_width = plot_width/2.54
plot_width = max(plot_width, 10)

if (nrow(comb) != 0) {
  hm_gene = hm_gene + baseMeanAnno_iso
  
  pdf("DE_Scaled_hm.pdf", width = plot_width, height = 12)
  draw(hm_gene)
  dev.off()
}