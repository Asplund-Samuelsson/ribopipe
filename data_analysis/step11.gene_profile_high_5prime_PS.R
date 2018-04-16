
### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input file
prof_file="results/2018-04-04/gene_PauseScore_profiles.with_plasmids.1.tab.gz"

### LOAD DATA ##################################################################

library(data.table)
gene_profiles = as.data.frame(fread(
  paste("gzip -dc ", prof_file, collapse=""),
  stringsAsFactors=F, header=T, sep="\t"
  ))

### DATA ANALYSIS ##############################################################

# Calculate distance from end of gene, and subset for 5 and 3 prime ends
gene_profiles$RelPosFromEnd = gene_profiles$RelPos - gene_profiles$Length + 1
gene_prof_5prime = subset(gene_profiles, RelPos %in% -50:-1)

# Subset to genes with high pause score
gene_prof_5prime_highPS = subset(gene_prof_5prime, PauseScore >= 1)

# How many genes?
high_PS_genes = unique(gene_prof_5prime_highPS$Name)

# It's 2474 genes, a bit much.

# Try different cutoffs
high_PS_genes = unique(subset(gene_prof_5prime, PauseScore >= 2)$Name) # 2058
high_PS_genes = unique(subset(gene_prof_5prime, PauseScore >= 1000)$Name) # 19

# Write those 18 genes to a file
write(high_PS_genes, "results/2018-04-04/high_5prime_PS_genes.with_plasmids.txt")

# Calculate average pause score in the 5' region
# If the pause score is close to one, it might indicate an extension
# of the coding/translated region
avg_PS_5p = aggregate(PauseScore ~ Sample + Name, gene_prof_5prime, mean)

# Try different ranges of average pause scores
candidate_genes = unique(subset(avg_PS_5p, PauseScore >= 0.99 & PauseScore <= 1.01)$Name)

# Write those 33 genes to a file
write(candidate_genes, "results/2018-04-04/similar_5prime_PS_genes.with_plasmids.txt")

# Look at the range up to -25 (that's where the RPM drops off)
avg_PS_5p_25 = aggregate(
  PauseScore ~ Sample + Name,
  subset(gene_prof_5prime, RelPos >= -25),
  mean
  )

candidate_genes = unique(subset(avg_PS_5p_25, PauseScore >= 10)$Name)

# Writing the 133 genes to a file
write(candidate_genes, "results/2018-04-04/25nt_5prime_PS_ge10_genes.with_plasmids.txt")
