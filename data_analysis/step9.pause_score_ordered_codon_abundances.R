#!/usr/bin/env Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
translation_file = args[1]
colour_file = args[2]

# Define infiles
codon_rpm_file = "analysis/codon_PauseScore_medians.tab"
codon_abundance_file = "analysis/count_codon.tab"

# Load data
codon = read.table(codon_rpm_file, header=T, sep="\t", stringsAsFactors=F)

t_table = read.table(translation_file, header=F, sep="\t", stringsAsFactors=F)
colnames(t_table) = c("Codon", "AminoAcid")

codon_abundance = read.table(codon_abundance_file, header=F, sep="\t", stringsAsFactors=F)
colnames(codon_abundance) = c("Count", "Codon")

# Add amino acids
codon = merge(codon, t_table)

# Add colours for amino acids
colours = scan(colour_file, character(), quote = "")
colours = c("#bdbdbd", colours, "#fc4e2a") # Add additional AA and STOP colours

amino_acids = c(grep("Ter", unique(codon$AminoAcid), value=T, invert=T), "Ter")

codon$AminoAcid = factor(codon$AminoAcid, levels=amino_acids)

# Plot it
library(ggplot2)
library(reshape2)

# Plot all samples with MADs
psMAD = codon

# Sort Codons by median pausescore
psMAD_median_means = aggregate(PauseScore ~ Codon, psMAD, mean)
psMAD_median_means = psMAD_median_means[order(psMAD_median_means$PauseScore),]
psMAD$Codon = factor(psMAD$Codon, levels=psMAD_median_means$Codon)

# Calculate max and min values
psMAD$Max = psMAD$PauseScore + psMAD$psMAD
psMAD$Min = psMAD$PauseScore - psMAD$psMAD

psMAD = merge(psMAD, t_table)
psMAD$AminoAcid = factor(psMAD$AminoAcid, levels=amino_acids)

# Sort Codons by median pausescore
psMAD_median_means = aggregate(PauseScore ~ Codon, psMAD, mean)
psMAD_median_means = psMAD_median_means[order(psMAD_median_means$PauseScore),]
psMAD$Codon = factor(psMAD$Codon, levels=psMAD_median_means$Codon)

# Add amino acids to abundances table
codon_abundance = merge(codon_abundance, t_table)

# Order Codons and AminoAcids the same way as in the psMAD table
codon_abundance$Codon = factor(codon_abundance$Codon, levels=levels(psMAD$Codon))
codon_abundance$AminoAcid = factor(codon_abundance$AminoAcid, levels=levels(psMAD$AminoAcid))

gp = ggplot(codon_abundance, aes(x=Codon, y=Count, fill=AminoAcid))
gp = gp + geom_col()
gp = gp + theme_bw()
gp = gp + scale_fill_manual(values=colours)
gp = gp + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7))

ggsave(
  "analysis/pause_score_ordered_codon_abundances.pdf",
  gp, width=240/25.4, height=100/25.4
  )
