#!/usr/bin/env Rscript

### LOAD DATA ##################################################################

# Define input file
args = commandArgs(trailingOnly=T)
indir = args[1]
annotation_file = args[2]
K = as.numeric(args[3])
# indir = "/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick/"
infile = normalizePath(file.path(indir, "analysis/gene_PauseScore_profiles.all_genes.tab.gz"))

# Load the file
library(data.table)
gene_profiles = as.data.frame(fread(
  paste("gzip -dc ", infile), sep="\t", stringsAsFactors=F
  ))

### ANALYZE ####################################################################

library(dplyr)

# Remove genes with fewer than 128 reads
gene_profiles = filter(gene_profiles, Reads >= 128)

# Clean up data
gene_profiles = filter(gene_profiles, is.finite(PauseScore))

# Subset to 5 prime UTR
gene_prof_5prime = filter(gene_profiles, RelPos %in% -50:-1)

# Drop uninteresting data
gene_prof_5prime = select(
  gene_prof_5prime, Name, Sample, Position, RPK_nucl, RelPos
)

# For neatness, use read counts instead
gene_prof_5prime = mutate(gene_prof_5prime, reads_nucl = RPK_nucl / 1000)

# Calculate RPK_UTR
gene_prof_5prime = group_by(gene_prof_5prime, Name, Sample) %>%
summarize(reads_UTR = sum(reads_nucl)) %>% as.data.frame() %>%
inner_join(gene_prof_5prime) %>% select(-RPK_nucl)

# Calculate RRO (Relative Ribosome Occupancy)
gene_prof_5prime = mutate(gene_prof_5prime, RRO = reads_nucl / reads_UTR)

# An RRO cannot be calculated for 5' UTRs with zero reads in total
gene_prof_5prime = filter(gene_prof_5prime, reads_UTR > 0)

### PLOT #######################################################################

library(ggplot2)

# Calculate averages
gene_prof_5prime_mean = group_by(gene_prof_5prime, Sample, RelPos) %>%
summarize(RRO = mean(RRO)) %>% as.data.frame()

# Mean lines of all genes, faceted by sample
gp = ggplot(gene_prof_5prime_mean, aes(x=RelPos, y=RRO))
gp = gp + geom_line(size=0.3)
gp = gp + theme_bw()
gp = gp + facet_grid(Sample ~ .)
gp = gp + scale_y_log10()
gp = gp + scale_x_continuous(minor_breaks=seq(-174, 168, 3))
gp = gp + theme(strip.background = element_blank())
gp = gp + annotation_logticks(sides="lr")

n_samples = length(unique(gene_prof_5prime_mean$Sample))

outfile = normalizePath(file.path(indir, "analysis/RRO_of_5prime_UTRs.pdf"))

ggsave(outfile, gp, height=6.875*n_samples*0.25, width=6.875, limitsize=F)


### CLUSTER PATTERNS ###########################################################

# Cast into a matrix
library(reshape2)

gene_prof_5prime_matrix = dcast(
  gene_prof_5prime, Name ~ RelPos + Sample, value.var="RRO"
)

rownames(gene_prof_5prime_matrix) = gene_prof_5prime_matrix$Name

gene_prof_5prime_matrix = select(gene_prof_5prime_matrix, -Name)

# Set missing values (UTR positions) to zero
gene_prof_5prime_matrix[is.na(gene_prof_5prime_matrix)] = 0

gene_prof_5prime_matrix = as.matrix(gene_prof_5prime_matrix)

library(pvclust)

# clst = pvclust(t(asinh(gene_prof_5prime_matrix)), parallel=T)
# Yields no significant clusters

library(bioDist)

clst = hclust(cor.dist(asinh(gene_prof_5prime_matrix)))
clst_cut = cutree(clst, k=K)

# Determine largest clusters
# clst_size = as.data.frame(table(clst_cut))
# clst_size = arrange(clst_size, -Freq)

# Take clusters > 20 genes
# clst_size_above_20 = filter(clst_size, Freq > 20)

# Extract genes per cluster
clst_cut_df = as.data.frame(clst_cut)
colnames(clst_cut_df) = "Cluster"
clst_cut_df$Name = rownames(clst_cut_df)

# clst_cut_df = filter(clst_cut_df, Cluster %in% clst_size_above_20$clst_cut)

# Plot those genes
gene_prof_5prime_clst = left_join(gene_prof_5prime, clst_cut_df)
#gene_prof_5prime_clst$Cluster[is.na(gene_prof_5prime_clst$Cluster)] = "Other"

gp = ggplot(gene_prof_5prime_clst, aes(x=RelPos, y=RRO, group=Name))
gp = gp + geom_line(size=0.3, alpha=0.2)
gp = gp + theme_bw()
gp = gp + facet_grid(Cluster ~ Sample)
#gp = gp + scale_y_log10()
gp = gp + scale_x_continuous(minor_breaks=seq(-174, 168, 3))
gp = gp + theme(strip.background = element_blank())
#gp = gp + annotation_logticks(sides="lr")

n_samples = length(unique(gene_prof_5prime_mean$Sample))

outfile = file.path(indir, "analysis/RRO_of_5prime_UTRs.clustered.pdf")

ggsave(outfile, gp, height=6.875*K*0.25, width=6.875*n_samples, limitsize=F)

# Find out what the genes in the clusters do...
# Load annotations
# annotation_file = "/home/johannes/proj/ribo/tools/ribopipe/genelists/syn_PCC6803/20180315_all_annotations_Micha.csv"
annotations = read.table(annotation_file, stringsAsFactors=F, sep=",", header=T)
colnames(annotations)[1] = "Name"

clst_process = left_join(gene_prof_5prime_clst, annotations) %>%
select(Name, Process, Cluster) %>% unique()
clst_process$Process[is.na(clst_process$Process)] = "No annotation"

# Add a cluster category with all clusters
clst_process = clst_process %>% mutate(Cluster = "All") %>% rbind(clst_process)

clst_process$Process = factor(
  clst_process$Process,
  levels=arrange(
    as.data.frame(table(filter(clst_process, Cluster == "All")$Process)),
    -Freq
  )$Var1
)

clst_process$Cluster = factor(
  clst_process$Cluster,
  levels=arrange(as.data.frame(table(clst_process$Cluster)), -Freq)$Var1
)

# Plot processes
gp = ggplot(clst_process, aes(x=Process))
gp = gp + geom_bar()
gp = gp + theme_bw()
gp = gp + facet_grid(~Cluster, scales="free_x")
gp = gp + coord_flip()
gp = gp + theme(strip.background = element_blank())

outfile = file.path(indir, "analysis/RRO_of_5prime_UTRs.processes_per_cluster.pdf")

ggsave(outfile, gp, height=6.875, width=6.875*(K+1)*0.25, limitsize=F)
