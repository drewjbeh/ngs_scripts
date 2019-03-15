library(ggplot2)
library(reshape2)
library(gridExtra)
library(plyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

mods = read.table("mods/mismatchTable.csv", header=T, sep = "\t", quote = '')
mods$proportion[is.na(mods$proportion)] = 0
mods$cluster = gsub(".*?_tRNA-", "", mods$cluster)
mods = subset(mods, bam == 'mm10_ESC/lib135_mt_E14_1_42_sr_on_trimFinal_SNP.unpaired_uniq.bam')
mods_agg = aggregate(mods$proportion, by = list(cluster=mods$cluster, pos=mods$pos, bam=mods$bam, struct=mods$struct), FUN = sum)
mods_wide = dcast(mods_agg[,c('cluster','pos', 'x')], list(.(cluster), .(pos)), value.var = 'x', fun.aggregate = mean)
mods_wide[is.na(mods_wide)] = 0
rownames(mods_wide) = mods_wide$cluster
mods_wide = mods_wide[, -1]
mods_mat = as.matrix(mods_wide)

stops = read.table("mods/RTstopTable.csv", header = T, sep = "\t", quote = '')
stops$proportion[is.na(stops$proportion)] = 0
stops$cluster = gsub(".*?_tRNA-", "", stops$cluster)
stops = subset(stops, bam == 'mm10_ESC/lib135_mt_E14_1_42_sr_on_trimFinal_SNP.unpaired_uniq.bam')
stops = stops[-5396, ] ## NB!! This cluster had info for pos 122 which doesn't exist - removed manually. Must be changed for other libraries ##
stops_wide = dcast(stops[,c('cluster','pos', 'proportion')], list(.(cluster), .(pos)), value.var = 'proportion', fun.aggregate = mean)
stops_wide[is.na(stops_wide)] = 0
rownames(stops_wide) = stops_wide$cluster
stops_wide = stops_wide[, -1]
stops_mat = as.matrix(stops_wide)

col_fun = colorRamp2(c(0, 0.5, 1), c("#f7fcf0", "#7bccc4", "#084081"))

col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(mods_agg$x, by = list(pos = mods_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
count_mods = mods_agg %>% group_by(cluster) %>% summarise(count = sum(x > 0.1))
row_anno = rowAnnotation(Count = row_anno_barplot(count_mods$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
mods_heatmap = Heatmap(mods_mat, row_title = "RT misincorporations", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "Misincorporation proportion", direction = "horizontal"))

col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(stops$proportion, by = list(pos = stops$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
count_stops = stops %>% group_by(cluster) %>% summarise(count = sum(proportion > 0.05))
row_anno = rowAnnotation(Count = row_anno_barplot(count_stops$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
stops_heatmap = Heatmap(stops_mat, row_title = "RT stops", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "RT stop proportion", direction = "horizontal"))

heatmap_list = mods_heatmap %v% stops_heatmap
draw(heatmap_list, ht_gap = unit(10, "mm"))

# Signatures of misinc (mouse)
# Select only mods of interest: m1G9, m2,2G26, m3C32, m1G37, m1A58
# ungapped alignment positions (seen from heatmaps or in mod context file): 13, 36, 42, 49, 91

context_info = read.table("mm10_modContext.txt", header = TRUE)
context_info$cluster = sub("Mus_musculus_tRNA-","", context_info$cluster)
sub_mods = subset(mods, pos %in% c(13, 36, 42, 49, 91))

# aggregate by type (sum) to get total misinc. proportions for each cluster
sub_mods_aggcluster = aggregate(sub_mods$proportion, by = list(cluster = sub_mods$cluster, pos = sub_mods$pos), FUN = sum)
cols = brewer.pal(9, "GnBu")[-(1:2)]
ggplot(sub_mods_aggcluster, aes(x=as.character(pos), y = x, color = x)) + geom_jitter(width = 0.1, size = 3) +
  theme_bw() + facet_wrap(~pos, scales = "free_x") + scale_color_gradientn(colours = cols) +
  geom_hline(yintercept = 0.1, linetype = "dashed", alpha = 0.4)

# create filter list of rows where total misinc. rate < 0.1 
filter_0.1 = subset(sub_mods_aggcluster, x < 0.1)

# filter these rows from mods table
sub_mods_aggtype = anti_join(sub_mods, filter_0.1, by=c("cluster","pos"))

# add in upstream nucl info
sub_mods_aggtype = merge(sub_mods_aggtype, context_info, by = c("cluster","pos"))

# aggregate accross all clusters to get mean misinc. pattern at each position
sub_mods_aggtype_up = aggregate(sub_mods_aggtype$proportion, by = list(pos = sub_mods_aggtype$pos, type = sub_mods_aggtype$type, upstream = sub_mods_aggtype$upstream), FUN = function(x) c(mean=mean(x), sd=sd(x)))
sub_mods_aggtype_up = do.call("data.frame", sub_mods_aggtype_up)
sub_mods_aggtype_down = aggregate(sub_mods_aggtype$proportion, by = list(pos = sub_mods_aggtype$pos, type = sub_mods_aggtype$type, downstream = sub_mods_aggtype$downstream), FUN = function(x) c(mean=mean(x), sd=sd(x)))
sub_mods_aggtype_down = do.call("data.frame", sub_mods_aggtype_down)
sub_mods_aggtype_all = aggregate(sub_mods_aggtype$proportion, by = list(identity = sub_mods_aggtype$identity, type = sub_mods_aggtype$type, upstream = sub_mods_aggtype$upstream, pos = sub_mods_aggtype$pos), FUN = function(x) c(mean=mean(x), sd=sd(x)))
sub_mods_aggtype_all = do.call("data.frame", sub_mods_aggtype_all)

# Signatures unstacked with error bars
# by identity, upstream and pos
ggplot(sub_mods_aggtype_all, aes(x = type, y = x.mean, fill = type)) + 
  geom_bar(stat="identity", width = 0.8, position =position_dodge(width=0.9), alpha = 0.9) + 
  facet_grid(upstream~pos+identity , scales = "free_x") + 
  geom_errorbar(aes(ymin = x.mean , ymax = x.mean + x.sd), width = 0.1, position = position_dodge(width=0.9)) + 
  theme_bw() +
  scale_fill_manual(values = c("#BB7A6C", "#7DB0A9", "#9F8FA9", "#8B9EB7"))

