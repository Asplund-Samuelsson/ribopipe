#!/usr/bin/env Rscript

### LOAD DATA ##################################################################

# Define input file
infile = "analysis/gene_PauseScore_profiles.tab.gz"
# infile = "/hdd/common/proj/RibosomeProfiling/results/2018-04-23/gene_PauseScore_profiles.NO_MULTI_ORF_POS.NO_ORF_BELOW_128R.tab.gz"

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
cl = makeCluster(16)
clusterExport(cl, c("lengths_of_5prime", "gene_prof_5prime", "rbindlist"))

sampled_5prime_PS = as.data.frame(rbindlist(parLapply(
  cl, 1:10000,
  function (n) {
    start = sample(lengths_of_5prime$MinRelPos, 1)
    end = -1
    random_5prime = as.data.frame(rbindlist(lapply(
      start:end,
      function (i) {
        aggregate(
          PauseScore ~ Sample,
          subset(gene_prof_5prime, RelPos == i),
          sample, 1
        )
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

cl = makeCluster(16)
clusterExport(cl, c("average_5prime_PS", "sampled_5prime_PS", "rbindlist"))

average_5prime_PS_test = as.data.frame(rbindlist(parLapply(
  cl, unique(average_5prime_PS$Name),
  function (gene) {
    high_random_5prime_PS = as.data.frame(rbindlist(lapply(
      c("A", "B", "C", "D", "E"),
      function (x) {
        subset(
          sampled_5prime_PS,
          Sample == x & PauseScore >= subset(
            average_5prime_PS, Sample == x & Name == gene
          )$PauseScore
        )
      }
    )))
    p_values = as.data.frame(table(high_random_5prime_PS$Sample) / 10000)
    colnames(p_values) = c("Sample", "p_val")
    merge(subset(average_5prime_PS, Name == gene), p_values)
  }
)))

stopCluster(cl)

# Adjust p-values per Sample
average_5prime_PS_test = as.data.frame(rbindlist(lapply(
  c("A", "B", "C", "D", "E"),
  function (x) {
    test_sub = subset(average_5prime_PS_test, Sample == x)
    test_sub$p_adj = p.adjust(test_sub$p_val, "BH")
    test_sub
  }
)))

# Adjusting the p-values leaves no significant genes

# Save list of significant genes
outfile = "analysis/significant_5prime_PS_genes.txt"
# outfile = "results/2018-04-23/significant_5prime_PS_genes.CSD2_seqmagick.2.txt"
write(
  unique(subset(average_5prime_PS_test, p_val < 0.05)$Name),
  outfile
)

# Save results table
outfile = "analysis/high_5prime_PS_test.tab"
# outfile = "results/2018-04-23/high_5prime_PS_test.CSD2_seqmagick.2.tab"
write.table(
  average_5prime_PS_test,
  outfile,
  quote=F, sep="\t", row.names=F, col.names=T
)
