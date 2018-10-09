#!/usr/bin/env Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
rp_file = args[1] # RP RPKM file
rn_file = args[2] # RNA RPKM file
# rp_file = "/hdd/common/proj/RibosomeProfiling/results/2018-08-30_CCE1_RP/2018.02.02.CCE1/analysis/all_CDS_RPKM_no_filter.tab"
# rn_file = "/hdd/common/proj/RNAseq/results/2018-09-14_CCE1/2018.06.11.CCE1_RNAseq/analysis/all_CDS_RPKM_no_filter.tab"

# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(scales)
library(tibble)

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
rp_wide = dcast(rp, Name ~ Sample, value.var="RPKM") %>% na.omit
rownames(rp_wide) = rp_wide$Name
rp_matrix = select(rp_wide, -Name) %>% as.matrix %>% t %>% log %>% scale %>% t

rn_wide = dcast(rn, Name ~ Sample, value.var="RPKM") %>% na.omit
rownames(rn_wide) = rn_wide$Name
rn_matrix = select(rn_wide, -Name) %>% as.matrix %>% t %>% log %>% scale %>% t

full_wide = rbind(rp, rn) %>% dcast(Name ~ Sample, value.var="RPKM") %>% na.omit
rownames(full_wide) = full_wide$Name
full_matrix = select(full_wide, -Name) %>% as.matrix %>% t %>% log %>% scale %>% t

# Perform PCA
rp_pca = prcomp(rp_matrix)
rn_pca = prcomp(rn_matrix)
full_pca = prcomp(full_matrix)

# Create plotting dataframes
plot_rp = as.data.frame(rp_pca$rotation)
plot_rp$Sample = rownames(plot_rp)

plot_rn = as.data.frame(rn_pca$rotation)
plot_rn$Sample = rownames(plot_rn)

plot_full = as.data.frame(full_pca$rotation)
plot_full$Sample = rownames(plot_full)

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

# Determine what genes have the greatest influence on the PCs
# library(factoextra)
# fviz_pca_var(
#   rp_pca,
#   col.var = "contrib", # Color by contributions to the PC
#   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#   repel = TRUE     # Avoid text overlapping
# )
#
# fviz_pca_biplot(
#   rp_pca, repel = TRUE,
#   col.var = "#2E9FDF", # Variables color
#   col.ind = "#696969"  # Individuals color
# )

# Extract genes in terms of PCs
rp_pca_x = as.data.frame(rp_pca$x) %>% select(PC1, PC2, PC3)
rp_pca_x$Name = rownames(rp_pca_x)

# Plot histograms of selected PCs
# hist(abs(rp_pca_x$PC1), 100)
# hist(abs(rp_pca_x$PC2), 100)
# hist(abs(rp_pca_x$PC3), 100)

# Calculate relative importance of genes for PCs 1 and 2
rp_pca_x = mutate(
  rp_pca_x,
  PC1_rel = abs(PC1) / sum(abs(PC1)), PC1_sign = sign(PC1),
  PC2_rel = abs(PC2) / sum(abs(PC2)), PC2_sign = sign(PC2)
)

# Determine if gene is above importance threshold for PCs 1 and 2
rp_pca_x = mutate(
  rp_pca_x,
  PC1_important = ifelse(PC1_rel > 0.00102, 1, 0),
  PC2_important = ifelse(PC2_rel > 0.00102, 1, 0)
)
# rp_pca_x = mutate(
#   rp_pca_x,
#   PC1_important = ifelse(abs(PC1) > 2.8, 1, 0),
#   PC2_important = ifelse(abs(PC2) > 2.3, 1, 0)
# )
# rp_pca_x = mutate(
#   rp_pca_x,
#   PC1_important = ifelse(PC1 < -1.3, 1, 0),
#   PC2_important = ifelse(PC2 < -1.2, 1, 0)
# )

# Create an Influence label for genes
rp_pca_x = mutate(
  rp_pca_x,
  Influence = ifelse(
    PC1_important + PC2_important == 2, "Both",
    ifelse(
      PC1_important & !PC2_important, "PC1",
      ifelse(
        !PC1_important & PC2_important, "PC2", "None"
      )
    )
  )
)

# Check if any genes influence both PCs
table(rp_pca_x$Influence)

# The genes with influence only have positive values for "their" PC

# Plot the expression patterns
rp_expr = inner_join(
  select(rp_pca_x, Name, Influence),
  mutate(as.data.frame(rp_matrix), Name = rownames(rp_matrix))
) %>% melt

rp_expr = inner_join(
  rp_expr,
  data.frame(
    variable = c(
      "A1_RPF", "A2_RPF", "B1_RPF", "B2_RPF",
      "C1_RPF", "C2_RPF", "D1_RPF", "D2_RPF",
      "E1_RPF"
    ),
    Time = c(
      -1, 23, 1, 25, 6, 30, 11, 35, 13
    )
  )
)

rp_expr = mutate(rp_expr, Z_score = value) %>% select(-value, variable)

# Add annotations
annotation_file = "/home/johannes/proj/ribo/tools/ribopipe/genelists/syn_PCC6803/20180315_all_annotations_Micha.csv"
annotations = read.table(annotation_file, stringsAsFactors=F, sep=",", header=T)
colnames(annotations)[1] = "Name"

rp_expr = left_join(rp_expr, annotations) %>%
  select(Name, Influence, Time, Z_score, Process)
rp_expr$Process[is.na(rp_expr$Process)] = "No annotation"

# Add relative PC contribution for colouring
rp_expr = left_join(
  rp_expr, rbind(
    rp_pca_x %>% filter(Influence == "PC1") %>%
      select(Name, Influence, PC1_rel) %>%
      mutate(PC_contribution = scale(PC1_rel)) %>% select(-PC1_rel),
    rp_pca_x %>% filter(Influence == "PC2") %>%
      select(Name, Influence, PC2_rel) %>%
      mutate(PC_contribution = scale(PC2_rel)) %>% select(-PC2_rel)
  )
)

gp = ggplot(rp_expr, aes(x=Time, y=Z_score, group=Name, colour=PC_contribution))
gp = gp + annotate(
  "rect", xmin=-1, xmax=0,
  ymin=min(rp_expr$Z_score)*1.02, ymax=max(rp_expr$Z_score)*1.02,
  fill="black", alpha=0.1
)
gp = gp + annotate(
  "rect", xmin=12, xmax=24,
  ymin=min(rp_expr$Z_score)*1.02, ymax=max(rp_expr$Z_score)*1.02,
  fill="black", alpha=0.1
)
gp = gp + annotate(
  "rect", xmin=36, xmax=37,
  ymin=min(rp_expr$Z_score)*1.02, ymax=max(rp_expr$Z_score)*1.02,
  fill="black", alpha=0.1
)
gp = gp + geom_line(alpha=0.2)
gp = gp + facet_grid(Process~Influence)
gp = gp + theme_bw()
gp = gp + theme(
  strip.background = element_blank(),
  strip.text.y = element_text(angle = 0, hjust=0)
)
gp = gp + scale_colour_gradient2(low="#542788", mid="#f7f7f7", high="#b35806")

outfile = paste(dirname(rp_file), "rp_PC_contribution_z_score_patterns.pdf", sep="/")
ggsave(outfile, gp, height=25/2.54, width=40/2.54)

# Save a table with the influence
rp_table_out = left_join(
  left_join(
    inner_join(
      unique(select(rp_expr, Name, Influence, Process)),
      rp %>% group_by(Name) %>% summarise(
        RPKM_mean = mean(RPKM),
        RPKM_SD = sd(RPKM),
        RPKM_CV = sd(RPKM)/mean(RPKM)*100
      ) %>% as.data.frame
    ),
    select(rp_pca_x, Name, PC1, PC2, PC1_rel, PC2_rel)
  ), annotations
)
outfile = paste(dirname(rp_file), "rp_PC_contribution_annotated.tab", sep="/")
write.table(rp_table_out, outfile, sep="\t", quote=F, row.names=F, col.names=T)

# Plot RPKM, CV, and contribution class
gp = ggplot(
  rp_table_out,
  aes(colour=log10(RPKM_mean), size=RPKM_CV, x=PC1, y=PC2, shape=Influence, alpha=Influence)
)
gp = gp + geom_vline(xintercept=min(filter(rp_table_out, Influence=="PC1")$PC1), colour="red")
gp = gp + geom_hline(yintercept=min(filter(rp_table_out, Influence=="PC2")$PC2), colour="red")
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + scale_alpha_discrete()
gp

################################################################################

# Use only the most varying genes
rp_scaled_long = melt(rp_matrix)
colnames(rp_scaled_long) = c("Name","Sample","Z_score")
rp_scaled_long_highz = filter(
  rp_scaled_long, Name %in% filter(rp_scaled_long, abs(Z_score) > 1.9)$Name
)
rp_wide_highz = dcast(rp_scaled_long_highz, Name ~ Sample)
rownames(rp_wide_highz) = rp_wide_highz$Name
rp_wide_highz = select(rp_wide_highz, -Name)
rp_matrix_highz = as.matrix(rp_wide_highz)

# Perform PCA
rp_pca_highz = prcomp(rp_matrix_highz)

# Create plotting dataframes
plot_rp_highz = as.data.frame(rp_pca_highz$rotation)
plot_rp_highz$Label = rownames(plot_rp_highz)

# Plot it
gp = ggplot(plot_rp_highz, aes(x=PC1, y=PC2, colour=PC3, label=Label))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + theme_bw()

gp

# Create plotting dataframes
plot_rp_highz = as.data.frame(rp_pca_highz$x)
plot_rp_highz$Name = rownames(plot_rp_highz)
plot_rp_highz = left_join(plot_rp_highz, annotations) %>%
  select(PC1, PC2, PC3, Name, Process)
plot_rp_highz$Process[is.na(plot_rp_highz$Process)] = "No annotation"

colour_file = "../main/art/2016-09-30/colours.txt"
colours = scan(colour_file, character(), quote = "")
# colours = c("#bdbdbd", colours, "#fc4e2a")

sample_arrows_rp_highz = as.data.frame(rp_pca_highz$rotation * 3)
sample_arrows_rp_highz$Sample = rownames(sample_arrows_rp_highz)


# Plot it
gp = ggplot(plot_rp_highz, aes(x=PC1, y=PC2, colour=Process, label=Name, group=Process))
gp = gp + geom_point()
gp = gp + scale_colour_manual(values=colours)
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + geom_segment(
  data=sample_arrows_rp_highz,
  mapping=aes(xend=PC1, yend=PC2, label=NA, group=NA),
  colour="#fc4e2a", x=0, y=0, arrow = arrow(length=unit(0.2, "cm"))
)
gp = gp + geom_text_repel(
  data=sample_arrows_rp_highz,
  mapping=aes(x=PC1, y=PC2, label=Sample, group=NA),
  colour="#fc4e2a", min.segment.length=10
)
#gp = gp + geom_density2d()
gp = gp + theme_bw()
gp = gp + coord_fixed()

outfile = paste(dirname(rp_file), "rp_PCA_clustering_of_genes.pdf", sep="/")

ggsave(outfile, gp, height=50/2.54, width=60/2.54)

################################################################################

# Add colours based on clusters from Jan
cluster_file = "/home/johannes/proj/ribo/data/2018-10-09/day-night_RPF_clusters_from_Jan.csv"
cl = read.table(cluster_file, sep=",", header=T, stringsAsFactors=F)
cl = select(cl, Name, Cluster = k4)
cl = mutate(cl, Cluster = as.character(Cluster))

plot_rp_highz = inner_join(plot_rp_highz, cl)

# Plot it
gp = ggplot(plot_rp_highz, aes(x=PC1, y=PC2, colour=Cluster, label=Name, group=Cluster))
gp = gp + geom_point()
gp = gp + scale_colour_manual(values=colours)
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + geom_segment(
  data=sample_arrows_rp_highz,
  mapping=aes(xend=PC1, yend=PC2, label=NA, group=NA),
  colour="#fc4e2a", x=0, y=0, arrow = arrow(length=unit(0.2, "cm"))
)
gp = gp + geom_text_repel(
  data=sample_arrows_rp_highz,
  mapping=aes(x=PC1, y=PC2, label=Sample, group=NA),
  colour="#fc4e2a", min.segment.length=10
)
#gp = gp + geom_density2d()
gp = gp + theme_bw()
gp = gp + coord_fixed()

outfile = paste(dirname(rp_file), "rp_PCA_clustering_of_genes.cluster_colors.pdf", sep="/")

ggsave(outfile, gp, height=40/2.54, width=48/2.54)

# Save one with all genes as well
rp_pca_x = as.data.frame(rp_pca$x) %>% select(PC1, PC2, PC3)
rp_pca_x$Name = rownames(rp_pca_x)

rp_pca_x = inner_join(rp_pca_x, cl)

sample_arrows_rp_pca_x = as.data.frame(rp_pca$rotation * 3)
sample_arrows_rp_pca_x$Sample = rownames(sample_arrows_rp_pca_x)


# Plot it
gp = ggplot(rp_pca_x, aes(x=PC1, y=PC2, colour=Cluster, label=Name, group=Cluster))
gp = gp + scale_colour_manual(values=colours)
gp = gp + geom_point(alpha=0.5)
gp = gp + geom_segment(
  data=sample_arrows_rp_pca_x,
  mapping=aes(xend=PC1, yend=PC2, label=NA, group=NA),
  colour="#fc4e2a", x=0, y=0, arrow = arrow(length=unit(0.2, "cm"))
)
gp = gp + geom_text_repel(
  data=sample_arrows_rp_pca_x,
  mapping=aes(x=PC1, y=PC2, label=Sample, group=NA),
  colour="#fc4e2a", min.segment.length=10
)
#gp = gp + geom_density2d()
gp = gp + theme_bw()
gp = gp + coord_fixed()

outfile = paste(dirname(rp_file), "rp_PCA_clustering_of_genes.cluster_colors.all_genes.pdf", sep="/")

ggsave(outfile, gp, height=25/2.54, width=30/2.54)
