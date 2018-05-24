#!/usr/bin/env Rscript

# Define infile
infile = "analysis/5prime_UTRs.min_6nt.right_ali.tab"

# Load table
utrs = read.table(infile, sep="\t", stringsAsFactors=F, header=F)
colnames(utrs) = c("Group", "Sequence")

# Split sequences into columns
library(data.table)

sequences = as.data.frame(rbindlist(lapply(
  strsplit(utrs$Sequence, ""), as.data.frame.list
)))

colnames(sequences) = paste("P", 0:49, sep="")

utrs = cbind(utrs[,"Group"], sequences)
colnames(utrs)[1] = "Group"

# Melt and recalculate position
library(reshape)
utrs = melt(utrs, id.vars="Group")
colnames(utrs)[2:3] = c("Position", "Nucleotide")
utrs$Position = as.numeric(sub("P", "", utrs$Position)) - 50

# Remove missing positions
utrs = subset(utrs, Nucleotide != "-")

# Make nucleotide factor
utrs$Nucleotide = factor(utrs$Nucleotide, levels = c("A","T","G","C"))

# Plot it
library(ggplot2)

gp = ggplot(
  utrs, aes(x=Position, fill=Nucleotide, group=Nucleotide)
)
gp = gp + geom_bar(position="dodge")
gp = gp + facet_grid(Group~Nucleotide, scales="free_y")
gp = gp + theme_bw()
gp = gp + scale_fill_manual(values=c("#80cdc1","#9970ab","#dfc27d","#2166ac"))
gp = gp + theme(strip.background = element_blank())

ggsave(
  "analysis/5prime_UTRs_nucleotide_count.pdf", gp,
  width = 360/25.4, height = 120/25.4
)
