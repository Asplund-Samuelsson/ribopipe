#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
RNAdir = args[1]

# Define input file
ribo_file="analysis/gene_RPKM.tab"

# List readsPerGene filenames
rna_files = list.files(
  path=paste(c(RNAdir, "/readsPerGene"), collapse=""),
  pattern="\\.readsPerGene\\.",
  full.names=T
  )

# Obtain sample and strand IDs for readsPerGene files in order
rna_samples = unlist(
  lapply(lapply(strsplit(rna_files, "\\."), tail, n=3L), head, n=1L)
  )
rna_strands = ifelse(
  unlist(lapply(strsplit(rna_files, "\\."), tail, n=1L)) == "p",
  "+",
  "-"
  )

# List totalNbrMappedReads filenames
rna_tc_files = list.files(
  path=paste(c(RNAdir, "/readcount"), collapse=""),
  pattern="totalNbrMappedReads",
  full.names=T
  )

# Obtain sample and strand IDs for readsPerGene files in order
rna_tc_samples = unlist(
  lapply(lapply(strsplit(rna_tc_files, "\\."), tail, n=2L), head, n=1L)
  )

### LOAD DATA ##################################################################

library(data.table)

# Load the files
rib = read.table(ribo_file, header=T, stringsAsFactors=F, sep="\t")
rna_data = lapply(rna_files, fread, header=F)
rna_tc_data = unlist(lapply(rna_tc_files, scan))

# Add columns with sample and strand IDs to each dataframe
for( i in seq_along(rna_data)){
  rna_data[[i]] = cbind(as.data.frame(rna_data[[i]]), Sample=rna_samples[i])
  rna_data[[i]] = cbind(as.data.frame(rna_data[[i]]), strand=rna_strands[i])
}

# Create one dataframe for rna data
rna = as.data.frame(rbindlist(rna_data))
colnames(rna)[1:5] = c("Sequence", "Name", "Start", "End", "Reads")

# Create total count data frames with sample names
tc_rna = data.frame(Sample = rna_tc_samples, TotalReads = rna_tc_data)

### NORMALIZE ##################################################################

# Add total count to per-gene RPM data
rna = merge(rna, tc_rna)

# Calculate RPM
rna$RPM = rna$Reads * 1000000 / rna$TotalReads

# Calculate gene length
rna$Length = rna$End - rna$Start + 1

# Special treatment for genes with negative start
rna$Length = ifelse(rna$Start > 0, rna$Length, rna$End - rna$Start)

# Calculate RPKM
rna$RPKM = 1000 * rna$RPM / rna$Length

# Combine RNAseq and ribosome profiling data
gen = rbind(
  rna[,grep("TotalReads", colnames(rna), invert=T)],
  rib[,grep("TotalReads", colnames(rib), invert=T)]
  )

# It is not interesting to include samples B, C, and D
gen = subset(gen, !(Sample %in% c("B","C","D")))

# Filter out low read count genes
# 3729 genes before
genes_below_threshold = unique(as.character(
  subset(gen, Reads < 32 & !(Sample %in% c("A","E")))$Name
  ))
gen = subset(gen, !(Name %in% genes_below_threshold))
# 2582 genes after

# Plot all vs all
library(GGally)
library(ggplot2)

library(reshape2)

gen$Name = gsub("'", '', gen$Name)
gen_wide = dcast(gen[,c("Sample","Name","RPKM")], Name ~ Sample, value.var="RPKM")

# Remove genes with no ribosome profiling data
gen_wide = subset(gen_wide, !is.na(A) & !is.na(E))

# Modify ggally_points function to add an abline
points_abline = function(data, mapping, ...){
  p = ggplot(data = data, mapping = mapping) +
      geom_point(..., alpha = 0.1, size=0.5) +
      geom_abline(slope=1, intercept=0, colour="red", alpha=0.6)
  p
}

# Plot all samples versus eachother
png("analysis/rnaseq_vs_ribprof.png", height = 800, width = 800)
g <- ggpairs(
  log10(gen_wide[,2:5]),
  lower = list(continuous = points_abline,
               combo = wrap("dot", alpha = 0.4, size=2)
               )
)
print(g)
dev.off()


# Check correlation of fold changes
gen_wide$fold_rib = gen_wide$E / gen_wide$A
gen_wide$fold_rna = gen_wide$time24 / gen_wide$time0

# Plot fold changes versus eachother
png("analysis/rnaseq_vs_ribprof.fold_change.png", height = 600, width = 600)
g <- ggpairs(
  log10(gen_wide[,6:7]),
  lower = list(continuous = points_abline,
               combo = wrap("dot", alpha = 0.4, size=2)
               )
)
print(g)
dev.off()


# Plot TE vs RNAseq
gen_wide$TE_t0 = gen_wide$A / gen_wide$time0
gen_wide$TE_t24 = gen_wide$E / gen_wide$time24

# Plot TE versus eachother
png("analysis/rnaseq_vs_ribprof.TE.png", height = 600, width = 600)
g <- ggpairs(
  log10(gen_wide[,8:9]),
  lower = list(continuous = points_abline,
               combo = wrap("dot", alpha = 0.4, size=2)
               )
)
print(g)
dev.off()

TE_data = melt(gen_wide[,c(1,8,9)])
colnames(TE_data) = c("Name", "Sample", "TE")
TE_data$Sample = ifelse(TE_data$Sample == "TE_t0", "0h", "24h")

# Plot TE distribution
gp = ggplot(TE_data, aes(x=TE, group=Sample, fill=Sample))
gp = gp + geom_density(alpha=0.3)
gp = gp + theme_bw()
gp = gp + scale_fill_manual(values=c("#1b7837","#762a83"))
gp = gp + scale_x_log10()
gp = gp + annotation_logticks(sides="b")

ggsave(
  "analysis/rnaseq_vs_ribprof.TE_distribution.pdf",
  gp, height = 90/25.4, width = 180/25.4
)

# Plot TE versus RNAseq
rna_plot = melt(gen_wide[,c(1,2,3)])
colnames(rna_plot) = c("Name", "Sample", "RNA_RPKM")
rna_plot$Sample = ifelse(rna_plot$Sample == "time0", "0h", "24h")

rna_plot = merge(rna_plot, TE_data)

gp = ggplot(rna_plot, aes(x=RNA_RPKM, y=TE, colour=Sample))
gp = gp + geom_point(size=0.4, alpha=0.3)
gp = gp + theme_bw()
gp = gp + facet_wrap(~Sample, ncol=2)
gp = gp + scale_x_log10()
gp = gp + scale_y_log10()
gp = gp + annotation_logticks(sides="bl")
gp = gp + scale_colour_manual(values=c("#1b7837","#762a83"))
gp = gp + theme(
  strip.background = element_blank()
)

ggsave(
  "analysis/rnaseq_vs_ribprof.TE_vs_RNAseq.pdf",
  gp, height = 90/25.4, width = 180/25.4
)

# What about fold change in TE vs fold change in RNAseq?
gen_wide$fold_TE = gen_wide$TE_t24 / gen_wide$TE_t0

# Plot fold changes versus eachother
png("analysis/rnaseq_vs_ribprof.TE_vs_RNAseq_fold_change.png", height = 600, width = 600)
g <- ggpairs(
  log10(gen_wide[,c(10,7)]),
  lower = list(continuous = points_abline,
               combo = wrap("dot", alpha = 0.4, size=2)
               )
)
print(g)
dev.off()

# Plot all ribosome profiling samples versus eachother
gen_wide = dcast(rib[,c("Sample","Name","RPKM")], Name ~ Sample, value.var="RPKM")

# Plot all samples versus eachother
png("analysis/ribo_all_vs_all_samples_and_genes.png", height = 800, width = 800)
g <- ggpairs(
  log(gen_wide[,2:ncol(gen_wide)]),
  lower = list(continuous = wrap("points", alpha = 0.1, size=0.5), combo = wrap("dot", alpha = 0.4, size=2) )
)
print(g)
dev.off()
