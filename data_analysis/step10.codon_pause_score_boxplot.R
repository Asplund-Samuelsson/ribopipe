#!/usr/bin/env Rscript

# Load command line arguments
args = commandArgs(trailingOnly=T)
translation_file = args[1]
colour_file = args[2]

# Define infiles
codon_rpm_file = "analysis/codon_PauseScore_medians.tab"

# Load data
codon = read.table(codon_rpm_file, header=T, sep="\t", stringsAsFactors=F)

t_table = read.table(translation_file, header=F, sep="\t", stringsAsFactors=F)
colnames(t_table) = c("Codon", "AminoAcid")

# Add amino acids
codon = merge(codon, t_table)

# Add colours for amino acids
colours = scan(colour_file, character(), quote = "")
colours = c("#bdbdbd", colours, "#fc4e2a") # Add additional AA and STOP colours

amino_acids = c(grep("Ter", unique(codon$AminoAcid), value=T, invert=T), "Ter")

codon$AminoAcid = factor(codon$AminoAcid, levels=amino_acids)

# Plot it
library(ggplot2)
library(ggrepel)

# Cast mean Pause Scores into wide format
library(reshape2)

# Sort Codons by median pausescore
codon_median_means = aggregate(PS_median ~ Codon, codon, mean)
codon_median_means = codon_median_means[order(codon_median_means$PS_median),]
codon$Codon = factor(codon$Codon, levels=codon_median_means$Codon)

# Add amino acids
codon = merge(codon, t_table)
codon$AminoAcid = factor(codon$AminoAcid, levels=amino_acids)

# asinh transformation and breaks for ggplot2
# http://wresch.github.io/2013/03/08/asinh-scales-in-ggplot2.html
library(scales)
asinh_breaks <- function(x) {
  br <- function(r) {
    lmin <- round(log10(r[1]))
    lmax <- round(log10(r[2]))
    lbreaks <- seq(lmin, lmax, by = 1)
    breaks <- 10 ^ lbreaks
  }
  p.rng <- range(x[x > 0], na.rm = TRUE)
  breaks <- br(p.rng)
  if (min(x) <= 0) {breaks <- c(0, breaks)}
  if (sum(x < 0) > 1) { #more negative values that expected from expanding scale that includes zero
    n.rng <- -range(x[x < 0], na.rm = TRUE)
    breaks <- c(breaks, -br(n.rng))
  }
  return(sort(breaks))
}
asinh_trans <- function() {
  trans_new("asinh",
            transform = asinh,
            inverse   = sinh,
            breaks = asinh_breaks)
}

custom_breaks = c(0, 10, 100)

gp = ggplot(codon,
    aes(
      x=Codon, ymin=PS_min, lower=PS_lower,
      middle=PS_median, upper=PS_upper,
      ymax=PS_max, colour=AminoAcid
    )
  )
gp = gp + geom_boxplot(stat="identity")
gp = gp + theme_bw()
gp = gp + facet_grid(Sample ~ .)
gp = gp + scale_colour_manual(values=colours)
gp = gp + scale_y_continuous(trans="asinh", breaks=custom_breaks, labels=custom_breaks)
#gp = gp + scale_y_sqrt()
gp = gp + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7))

ggsave(
  "analysis/codon_pause_score_boxplot.pdf",
  gp, width=240/25.4, height=180/25.4
  )
