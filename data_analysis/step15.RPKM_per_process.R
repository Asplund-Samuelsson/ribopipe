#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
genelist_file = args[1]
annotation_file = args[2]
indir = "." # A ribopipe results directory
# genelist_file = "/home/johannes/proj/ribo/tools/ribopipe/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"
# annotation_file = "/home/johannes/proj/ribo/tools/ribopipe/genelists/syn_PCC6803/20180315_all_annotations_Micha.csv"

# List RPM0 filenames
rpm_files = list.files(
  path=paste(c(indir, "/RPM"), collapse=""), pattern="\\.RPM0\\.", full.names=T
  )

# List readsPerGene filenames
gene_files = list.files(
  path=paste(c(indir, "/readsPerGene"), collapse=""),
  pattern="\\.readsPerGene\\.",
  full.names=T
  )

# Obtain sample and strand IDs for RPM0 files in order
rpm_samples = unlist(
  lapply(lapply(strsplit(rpm_files, "\\."), tail, n=3L), head, n=1L)
  )
rpm_strands = ifelse(
  unlist(lapply(strsplit(rpm_files, "\\."), tail, n=1L)) == "p",
  "+",
  "-"
  )

# Obtain sample and strand IDs for readsPerGene files in order
gene_samples = unlist(
  lapply(lapply(strsplit(gene_files, "\\."), tail, n=3L), head, n=1L)
  )
gene_strands = ifelse(
  unlist(lapply(strsplit(gene_files, "\\."), tail, n=1L)) == "p",
  "+",
  "-"
  )

# List totalNbrMappedReads filenames
gen_tc_files = list.files(
  path=paste(c(indir, "/readcount"), collapse=""),
  pattern="totalNbrMappedReads",
  full.names=T
  )

# Obtain sample and strand IDs for readsPerGene files in order
gen_tc_samples = unlist(
  lapply(lapply(strsplit(gen_tc_files, "\\."), tail, n=2L), head, n=1L)
  )

### LOAD DATA ##################################################################

library(data.table)

# Load the files
rpm_data = lapply(rpm_files, fread, header=F)
gene_data = lapply(gene_files, fread, header=F)
gen_tc_data = unlist(lapply(gen_tc_files, scan))

# Add columns with sample and strand IDs to each dataframe
for( i in seq_along(rpm_data)){
  rpm_data[[i]] = cbind(as.data.frame(rpm_data[[i]]), Sample=rpm_samples[i])
  rpm_data[[i]] = cbind(as.data.frame(rpm_data[[i]]), strand=rpm_strands[i])
}
for( i in seq_along(gene_data)){
  gene_data[[i]] = cbind(as.data.frame(gene_data[[i]]), Sample=gene_samples[i])
  gene_data[[i]] = cbind(as.data.frame(gene_data[[i]]), strand=gene_strands[i])
}

# Create one dataframe for RPM data
rpm = as.data.frame(rbindlist(rpm_data))
colnames(rpm)[1:3] = c("Sequence", "Position", "RPM")

# Create one dataframe for gene data
gen = as.data.frame(rbindlist(gene_data))
colnames(gen)[1:5] = c("Sequence", "Name", "Start", "End", "Reads")

# Create total count data frames with sample names
tc_gen = data.frame(Sample = gen_tc_samples, TotalReads = gen_tc_data)

# Load gene types
gene_types = read.table(
  genelist_file, stringsAsFactors=F, sep="\t", header=T, quote=""
)[,c("Type","Old_locus_tag", "Sequence")]
colnames(gene_types)[grep("Old", colnames(gene_types))] = "Name"

# Load annotations
annotations = read.table(annotation_file, stringsAsFactors=F, sep=",", header=T)
colnames(annotations)[1] = "Name"

### SHIFT RPM VALUES ###########################################################

# SIGNAL SHIFT
x = -12

# For the plus strand, the position shift is added
# For the minus strand, the position shift is subtracted

# Change position by the shift
rpm_shift = rpm
rpm_shift$Shift = x
rpm_shift$Position = ifelse(
  rpm_shift$strand == "+",
  rpm_shift$Position + rpm_shift$Shift,
  rpm_shift$Position - rpm_shift$Shift
  )

# Fold around beginning of circular genome
library(dplyr)

sequence_sizes = aggregate(Position ~ Sequence, rpm, max)
colnames(sequence_sizes)[2] = "Size"

rpm_shift = inner_join(rpm_shift, sequence_sizes)

rpm_shift$Position = ifelse(
  rpm_shift$Position < 1,
  rpm_shift$Position + rpm_shift$Size,
  ifelse(
    rpm_shift$Position > rpm_shift$Size,
    rpm_shift$Position - rpm_shift$Size,
    rpm_shift$Position
    )
  )

### RECALCULATE RPM PER GENE ###################################################

# Remove current read count
gen = gen[,grep("(Sample)|(Reads)", colnames(gen), invert=T, perl=T)]
gen = unique(gen)

# Subset to coding sequences only
gen = subset(gen, Name %in% subset(gene_types, Type == "CDS")$Name)

# Expand each gene to all positions
gen_allpos = as.data.frame(rbindlist(lapply(
  gen$Name,
  function (x) {
    gen_s = subset(gen, Name == x)
    Position = (gen_s$Start):(gen_s$End)
    merge(gen_s, data.frame(Position = Position))
    }
)))

# Fold around beginning of circular genome
gen_allpos = inner_join(gen_allpos, sequence_sizes)

gen_allpos$Position = ifelse(
  gen_allpos$Position < 1,
  gen_allpos$Position + gen_allpos$Size,
  ifelse(
    gen_allpos$Position > gen_allpos$Size,
    gen_allpos$Position - gen_allpos$Size,
    gen_allpos$Position
    )
  )

# Add shifted RPM data
gen_allpos = inner_join(gen_allpos, rpm_shift)

### NORMALIZE ##################################################################

# Add total reads
gen_allpos = inner_join(gen_allpos, tc_gen)

# Recalculate reads per nucleotide
gen_allpos$Reads = gen_allpos$RPM * gen_allpos$TotalReads / 1000000

# Calculate total reads on ORFs
total_ORF_reads = aggregate(Reads ~ Sample, gen_allpos, sum)
colnames(total_ORF_reads)[2] = "TotalReadsORF"
gen_allpos = inner_join(gen_allpos, total_ORF_reads)

# Recalculate RPM
gen_allpos$RPM = gen_allpos$Reads / gen_allpos$TotalReadsORF * 1000000

# Calculate the total reads of the genes
gen = inner_join(gen, aggregate(Reads ~ Name + Sample + TotalReadsORF, gen_allpos, sum))

# Calculate the total RPM of the genes
gen$RPM = gen$Reads / gen$TotalReadsORF * 1000000

# Remove entries that have less than 128 reads for each sample
gen = subset(gen, Reads > 128)

# Calculate gene length
gen$Length = gen$End - gen$Start + 1

# Special treatment for genes with negative start
gen$Length = ifelse(gen$Start > 0, gen$Length, gen$End - gen$Start)

# Calculate RPKM
gen$RPKM = 1000 * gen$RPM / gen$Length


### PLOT PROCESS BARS ##########################################################

# Merge
gen = merge(gen, annotations[,c("Name","Process")], all.x = T)

gen$Process[is.na(gen$Process)] = "No annotation"

# Clean up
gen = gen[,c("Name","Sample","RPKM","Process")]

# Sum per process
gen_sum = aggregate(RPKM ~ Sample + Process, gen, sum)
gen_sum$Sample = factor(gen_sum$Sample, levels = rev(c("A","B","C","D","E")))

# Sort processes by total RPKM of genes
sum_rpkm = aggregate(RPKM ~ Process, gen, sum)
sum_rpkm = sum_rpkm[order(sum_rpkm[,2]),]
gen_sum$Process = factor(
  gen_sum$Process,
  levels = sum_rpkm$Process
)

# Plot it
library(ggplot2)

gp = ggplot(gen_sum, aes(x=Process, colour=Sample, y=RPKM))
gp = gp + geom_point(alpha=0.6)
gp = gp + scale_colour_manual(values=rev(c("#1b7837","#5aae61","#a6dba0","#d9f0d3","#762a83")))
gp = gp + theme_bw()
gp = gp + coord_flip()
gp = gp + theme(
  legend.position = c(-1.15, 0.75),
  strip.background = element_blank()
)
gp = gp + guides(colour = guide_legend(reverse = T))

ggsave(
  "analysis/RPKM_per_process_points.pdf", gp,
  width = 180/25.4, height = 85/25.4
)

# Add information about 5' UTR
significant_file = "analysis/high_5prime_PS_test.tab"

average_5prime_PS = read.table(
  significant_file, header=T, sep="\t",stringsAsFactors=F
)

average_5prime_PS = mutate(
  average_5prime_PS,
  UTR = ifelse(PauseScore > PS_0.05, "High PS 5' UTR", "Other 5' UTR")
)

gen_utr = inner_join(gen, average_5prime_PS[,c("Name","Sample","UTR")])

# Sum per process
gen_sum = aggregate(RPKM ~ Sample + Process + UTR, gen_utr, sum)
gen_sum$Sample = factor(gen_sum$Sample, levels = rev(c("A","B","C","D","E")))

gen_sum$Process = factor(
  gen_sum$Process,
  levels = sum_rpkm$Process
)

# Plot
gp = ggplot(gen_sum, aes(x=Process, colour=Sample, y=RPKM))
gp = gp + geom_point(alpha=0.6)
gp = gp + scale_colour_manual(values=rev(c("#1b7837","#5aae61","#a6dba0","#d9f0d3","#762a83")))
gp = gp + theme_bw()
gp = gp + coord_flip()
gp = gp + facet_wrap(~UTR, scales="free_x")
gp = gp + theme(
  legend.position = c(-0.50, 0.8),,
  strip.background = element_blank()
)
gp = gp + guides(colour = guide_legend(reverse = T))

ggsave(
  "analysis/RPKM_per_process_points_by_UTR.pdf", gp,
  width = 280/25.4, height = 85/25.4
)
