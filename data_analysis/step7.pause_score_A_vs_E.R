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

# Adding regression line
# https://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph
lm_eqn <- function(df){
  m <- lm(A ~ E, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
     list(a = format(coef(m)[1], digits = 2),
          b = format(coef(m)[2], digits = 2),
         r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Cast mean Pause Scores into wide format
library(reshape2)
pause_score_wide = dcast(
  codon[,c("Codon","Sample","PauseScore")],
  Codon ~ Sample, value.var = "PauseScore"
  )

# Add amino acids
pause_score_wide = merge(pause_score_wide, t_table)
pause_score_wide$AminoAcid = factor(pause_score_wide$AminoAcid, levels=amino_acids)

gp = ggplot(pause_score_wide, aes(x=E, y=A, colour=AminoAcid))
gp = gp + geom_smooth(method = "lm", color="black", formula = y ~ x)
gp = gp + geom_point()
gp = gp + geom_text_repel(mapping=aes(label=Codon), size=3)
gp = gp + geom_text(x = 0.15, y = 0.5, label = lm_eqn(pause_score_wide), colour="black", parse = TRUE)
gp = gp + theme_bw()
gp = gp + scale_colour_manual(values=colours)
gp = gp + ggtitle("A (0h) vs E (24h) Median Codon Pause Score (-12 shift)")

ggsave("analysis/pause_score_A_vs_E.codon_medians.pdf", gp, width=180/25.4, height=120/25.4)


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

gp = ggplot(psMAD, aes(x=Codon, y=PauseScore, ymax=Max, ymin=Min, colour=AminoAcid))
gp = gp + geom_errorbar()
gp = gp + geom_line()
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + facet_grid(Sample ~ .)
gp = gp + scale_colour_manual(values=colours)
gp = gp + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7))

ggsave(
  "analysis/pause_score_distribution.codon_medians.pdf",
  gp, width=240/25.4, height=180/25.4
  )
