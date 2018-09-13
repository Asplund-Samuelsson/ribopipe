#!/usr/bin/env Rscript

# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(scales)

# Load command line arguments
args = commandArgs(trailingOnly=T)
rp_file = args[1] # RP RPKM file
rn_file = args[2] # RNA RPKM file
# rp_file = "/hdd/common/proj/RibosomeProfiling/results/2018-08-30_CCE1_RP/2018.02.02.CCE1/analysis/all_CDS_RPKM_no_filter.tab"
# rn_file = "/ssd/jan/seqdata/2018.06.11.CCE1_RNAseq/analysis/all_CDS_RPKM_no_filter.tab"

# Load data
rp = read.table(rp_file, sep="\t", stringsAsFactors=F, header=T) %>%
  mutate(Source = "RPF") %>%
  mutate(Sample = paste(Sample, Source, sep="_")) %>%
  select(Name, Sample, Source, Reads, RPKM) %>%
  filter(Reads >= 128)

rn = read.table(rn_file, sep="\t", stringsAsFactors=F, header=T) %>%
  mutate(Source = "RNA") %>%
  mutate(Sample = paste(Sample, Source, sep="_")) %>%
  select(Name, Sample, Source, Reads, RPKM) %>%
  filter(Reads >= 32)

# Log-transform and scale the data
rp_matrix = dcast(rp, Name ~ Sample, value.var="RPKM") %>%
  na.omit %>% select(-Name) %>% as.matrix %>% t %>% log %>% scale %>% t

rn_matrix = dcast(rn, Name ~ Sample, value.var="RPKM") %>%
  na.omit %>% select(-Name) %>% as.matrix %>% t %>% log %>% scale %>% t

# Perform PCA
rp_pca = prcomp(rp_matrix)
rn_pca = prcomp(rn_matrix)

# Create plotting dataframes
plot_rp = as.data.frame(rp_pca$rotation)
plot_rp$Sample = rownames(plot_rp)

plot_rn = as.data.frame(rn_pca$rotation)
plot_rn$Sample = rownames(plot_rn)

# Plot it

gp = ggplot(plot_rp, aes(x=PC1, y=PC2, colour=PC3, label=Sample))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + theme_bw()

outfile = paste(dirname(rp_file), "rp_3D-PCA_clustering_of_samples.pdf", sep="/")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)

# Calculate fraction of variance per PC
message("RPF PCs variance explained")
percent(rp_pca$sdev^2 / sum(rp_pca$sdev^2))

gp = ggplot(plot_rn, aes(x=PC1, y=PC2, colour=PC3, label=Sample))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + theme_bw()

outfile = paste(dirname(rp_file), "rn_3D-PCA_clustering_of_samples.pdf", sep="/")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)

# Calculate fraction of variance per PC
message("RNA PCs variance explained")
percent(rn_pca$sdev^2 / sum(rn_pca$sdev^2))


# Now perform PCA with all samples together
full_matrix = rbind(rp, rn) %>% dcast(Name ~ Sample, value.var="RPKM") %>%
  na.omit %>% select(-Name) %>% as.matrix %>% t %>% log %>% scale %>% t

# Perform PCA
full_pca = prcomp(full_matrix)

# Create plotting dataframes
plot_full = as.data.frame(full_pca$rotation)
plot_full$Sample = rownames(plot_full)

# Plot it
gp = ggplot(plot_full, aes(x=PC1, y=PC2, colour=PC3, label=Sample))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + theme_bw()

outfile = paste(dirname(rp_file), "rp+rn_3D-PCA_clustering_of_samples.pdf", sep="/")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)

# Calculate fraction of variance per PC
message("RP+RNA PCs variance explained")
percent(full_pca$sdev^2 / sum(full_pca$sdev^2))
