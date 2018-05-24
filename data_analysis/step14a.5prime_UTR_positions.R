#!/usr/bin/env Rscript

# Define input file
infile = "analysis/gene_PauseScore_profiles.tab.gz"
significant_file = "analysis/significant_5prime_PS_genes.txt"
# infile = "/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick/analysis/gene_PauseScore_profiles.tab.gz"
# significant_file = "/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick/analysis/significant_5prime_PS_genes.txt"

# Load the file
library(data.table)
library(dplyr)
gene_profiles = as.data.frame(fread(
  paste("gzip -dc ", infile), sep="\t", stringsAsFactors=F
  ))
significant_genes = scan(significant_file, character())

# Subset for 5 prime end
gene_prof_5prime = filter(
  gene_profiles,
  RelPos %in% -50:-1 & Name %in% significant_genes
)
gene_prof_5prime = unique(gene_prof_5prime[,
  c( "Name", "Sequence", "strand", "Position")
])

# Write to file
outfile = "analysis/significant_5prime_UTR_positions.tab"
write.table(
  gene_prof_5prime,
  outfile,
  quote=F, sep="\t", row.names=F, col.names=T
)

# Subset for 5 prime end
gene_prof_5prime = filter(
  gene_profiles,
  RelPos %in% -50:-1 & !(Name %in% significant_genes)
)
gene_prof_5prime = unique(gene_prof_5prime[,
  c( "Name", "Sequence", "strand", "Position")
])

# Write to file
outfile = "analysis/insignificant_5prime_UTR_positions.tab"
write.table(
  gene_prof_5prime,
  outfile,
  quote=F, sep="\t", row.names=F, col.names=T
)
