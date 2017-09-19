
### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input directory
indir="/ssd/common/proj/RibosomeProfiling/tsspipe_testing/2017.09.05.CSD2/readsPerGene"

# List readsPerGene filenames
gene_files = list.files(path=indir, pattern="\\.readsPerGene\\.", full.names=T)

# Obtain sample and strand IDs for readsPerGene files in order
gene_samples = unlist(
  lapply(lapply(strsplit(gene_files, "\\."), tail, n=3L), head, n=1L)
  )
gene_strands = ifelse(
  unlist(lapply(strsplit(gene_files, "\\."), tail, n=1L)) == "p",
  "+",
  "-"
  )

### LOAD DATA ##################################################################

# Load the files
gene_data = lapply(gene_files, read.delim, header=F)

# Add columns with sample and strand IDs to each dataframe
for( i in seq_along(gene_data)){
  gene_data[[i]] = cbind(gene_data[[i]], Sample=gene_samples[i])
  gene_data[[i]] = cbind(gene_data[[i]], Strand=gene_strands[i])
}

# Create one dataframe for gene data
gen = do.call("rbind", gene_data)
colnames(gen)[1:4] = c("Name", "Start", "End", "Reads")

# Calculate gene length
gen$Length = gen$End - gen$Start + 1

# Count total reads mapped to genes
total_count = aggregate(Reads ~ Sample, gen, sum)
colnames(total_count)[2] = "TotalReads"

# Merge with gen to create new column with total read count per sample
gen = merge(gen, total_count)

# Calculate RPKM
gen$RPKM = gen$Reads * 1000 * 1000000 / gen$Length / gen$TotalReads

# Filter out low read count genes
# 3317 genes before
genes_below_threshold = as.character(unique(subset(gen, Reads <= 128)$Name))
gen = subset(gen, !(Name %in% genes_below_threshold))

# Plot all vs all
library(GGally)
library(ggplot2)

library(reshape2)

gen$Name = gsub("'", '', gen$Name)
gen_wide = dcast(gen[,c("Sample","Name","RPKM")], Name ~ Sample, value.var="RPKM")

# Save Wide and Long format of RPKM values
write.table(gen_wide, "/tmp/CSD2.RPKM.wide.tab", sep="\t", quote=F, row.names=F, col.names=T)
write.table(gen, "/tmp/CSD2.RPKM.long.tab", sep="\t", quote=F, row.names=F, col.names=T)

# Plot all samples versus eachother
png("/tmp/ribo_all_vs_all_samples_and_genes.png", height = 800, width = 800)
g <- ggpairs(
  log(gen_wide[,2:6]),
  lower = list(continuous = wrap("points", alpha = 0.1, size=0.5), combo = wrap("dot", alpha = 0.4, size=2) )
)
print(g)
dev.off()
