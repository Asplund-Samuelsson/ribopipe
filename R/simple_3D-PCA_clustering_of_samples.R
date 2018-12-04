#!/usr/bin/env Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
rpkm_file = args[1] # RPKM file
# rpkm_file = "/hdd/common/proj/RNAseq/results/2018-09-14_CCE1/2018.06.11.CCE1_RNAseq/analysis/all_CDS_RPKM_no_filter.tab"

# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(scales)
library(tibble)

# Load data
rn = read.table(rpkm_file, sep="\t", stringsAsFactors=F, header=T) %>%
  select(Name, Sample, Reads, RPKM) %>%
  filter(Reads >= 32)

# Log-transform and scale the data
rn_wide = dcast(rn, Name ~ Sample, value.var="RPKM") %>% na.omit
rownames(rn_wide) = rn_wide$Name
rn_matrix = select(rn_wide, -Name) %>% as.matrix %>% t %>% log %>% scale %>% t

# Perform PCA
rn_pca = prcomp(rn_matrix)

# Create plotting dataframes
plot_rn = as.data.frame(rn_pca$rotation)
plot_rn$Sample = rownames(plot_rn)

# Calculate fraction of variance per PC
var_pc = percent(rn_pca$sdev^2 / sum(rn_pca$sdev^2))[1:3]

gp = ggplot(plot_rn, aes(x=PC1, y=PC2, colour=PC3, label=Sample))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + labs(
            x=paste("PC1 (", var_pc[1], ")", sep=""),
            y=paste("PC2 (", var_pc[2], ")", sep=""),
            colour=paste("PC3 (", var_pc[3], ")", sep="")
          )
gp = gp + theme_bw()

outfile = paste(dirname(rpkm_file), "3D-PCA_clustering_of_samples.pdf", sep="/")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)
