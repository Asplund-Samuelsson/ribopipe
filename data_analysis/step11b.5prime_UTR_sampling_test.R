#!/usr/bin/env Rscript

### LOAD DATA ##################################################################

# Define input file
infile = "analysis/gene_PauseScore_profiles.tab.gz"
# infile = "/hdd/common/proj/RibosomeProfiling/results/2018-04-23/gene_PauseScore_profiles.NO_MULTI_ORF_POS.NO_ORF_BELOW_128R.tab.gz"
# infile = "/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick/analysis/gene_PauseScore_profiles.tab.gz"

# Load the file
library(data.table)
library(dplyr)
gene_profiles = as.data.frame(fread(
  paste("gzip -dc ", infile), sep="\t", stringsAsFactors=F
  ))

# Subset for 5 prime end
gene_prof_5prime = filter(gene_profiles, RelPos %in% -50:-1)
gene_prof_5prime = gene_prof_5prime[,
  c("Name", "Sample", "RelPos", "PauseScore")
]

### TEST #######################################################################

# Calculate the distribution of 5' UTR lengths
lengths_of_5prime = aggregate(
  RelPos~Name, filter(gene_prof_5prime, Sample == "A"), min
  )
colnames(lengths_of_5prime)[2] = "MinRelPos"

# All ORFs with a minimum RelPos >= 0 have been removed above

# Calculate the average PS of each real UTR
average_5prime_PS = aggregate(
  PauseScore ~ Sample + Name, gene_prof_5prime, mean
  )

# Create list of random 5' UTRs
library(parallel)
cl = makeCluster(32)
clusterExport(
  cl,
  c("lengths_of_5prime", "gene_prof_5prime", "rbindlist", "filter", "data.table")
  )

sampled_5prime_PS = as.data.frame(rbindlist(parLapply(
  cl, 1:10000,
  function (n) {
    start = sample(lengths_of_5prime$MinRelPos, 1)
    end = -1
    random_5prime = as.data.frame(rbindlist(lapply(
      start:end,
      function (i) {
        as.data.frame(data.table(filter(gene_prof_5prime, RelPos == i))[,
          list(PauseScore = sample(PauseScore, 1)), by = 'Sample'
        ])
      }
    )))
    merge(
      aggregate(PauseScore ~ Sample, random_5prime, mean),
      data.frame(Iteration = n)
    )
  }
)))

stopCluster(cl)

# Compare the real UTR average PS to the random list
p_limits = aggregate(
  PauseScore ~ Sample, sampled_5prime_PS, quantile, c(0.95, 0.99, 0.999)
)

p_limits$PS_0.05 = p_limits$PauseScore[,"95%"]
p_limits$PS_0.01 = p_limits$PauseScore[,"99%"]
p_limits$PS_0.001 = p_limits$PauseScore[,"99.9%"]
p_limits = p_limits[,c(1,3:5)]

# Add p limits to average 5prime PS
average_5prime_PS = merge(average_5prime_PS, p_limits[,c("Sample","PS_0.05")])

# Save list of significant genes
outfile = "analysis/significant_5prime_PS_genes.txt"
write(
  unique(subset(average_5prime_PS, PauseScore > PS_0.05)$Name),
  outfile
)

# Save results table
outfile = "analysis/high_5prime_PS_test.tab"
write.table(
  average_5prime_PS,
  outfile,
  quote=F, sep="\t", row.names=F, col.names=T
)
