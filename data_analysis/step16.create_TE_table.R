#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Load command line arguments
args = commandArgs(trailingOnly=T)

# Define RNA directory
RNAdir = args[1]
# RNAdir = "/hdd/common/proj/RNAseq/results/2018-03-29/CO2_lim_with_plasmids"

# Define gene RPKM file
ribo_file = "analysis/gene_RPKM.tab"
# ribo_file="/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick/analysis/gene_RPKM.tab"

# Define annotations file
genelist_file = args[2]
annotation_file = args[3]
# genelist_file = "tools/ribopipe/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"
# annotation_file = "data/2018-05-03/20180315_all_annotations_Micha.csv"

# List readsPerGene filenames
rna_files = list.files(
  path=paste(c(RNAdir, "/readsPerGene"), collapse=""),
  pattern="\\.readsPerGene\\.",
  full.names=T
  )

# Obtain sample and strand IDs for readsPerGene files in order
rna_samples = unlist(
  lapply(lapply(strsplit(rna_files, "\\."), tail, n=3L), head, n=1L)
  )
rna_strands = ifelse(
  unlist(lapply(strsplit(rna_files, "\\."), tail, n=1L)) == "p",
  "+",
  "-"
  )

# List totalNbrMappedReads filenames
rna_tc_files = list.files(
  path=paste(c(RNAdir, "/readcount"), collapse=""),
  pattern="totalNbrMappedReads",
  full.names=T
  )

# Obtain sample and strand IDs for readsPerGene files in order
rna_tc_samples = unlist(
  lapply(lapply(strsplit(rna_tc_files, "\\."), tail, n=2L), head, n=1L)
  )

### LOAD DATA ##################################################################

library(data.table)

# Load the files
rib = read.table(ribo_file, header=T, stringsAsFactors=F, sep="\t")
rna_data = lapply(rna_files, fread, header=F)
rna_tc_data = unlist(lapply(rna_tc_files, scan))

# Add columns with sample and strand IDs to each dataframe
for( i in seq_along(rna_data)){
  rna_data[[i]] = cbind(as.data.frame(rna_data[[i]]), Sample=rna_samples[i])
  rna_data[[i]] = cbind(as.data.frame(rna_data[[i]]), strand=rna_strands[i])
}

# Create one dataframe for rna data
rna = as.data.frame(rbindlist(rna_data))
colnames(rna)[1:5] = c("Sequence", "Name", "Start", "End", "Reads")

# Create total count data frames with sample names
tc_rna = data.frame(Sample = rna_tc_samples, TotalReads = rna_tc_data)

# Load annotations
annotations = read.table(annotation_file, stringsAsFactors=F, sep=",", header=T)
colnames(annotations)[1] = "Name"
gene_types = read.table(
  genelist_file, stringsAsFactors=F, sep="\t", header=T, quote=""
)[,c("Type","Old_locus_tag", "Sequence")]
colnames(gene_types)[grep("Old", colnames(gene_types))] = "Name"

### NORMALIZE ##################################################################

# Add total count to per-gene RPM data
rna = merge(rna, tc_rna)

# Calculate RPM
rna$RPM = rna$Reads * 1000000 / rna$TotalReads

# Calculate gene length
rna$Length = rna$End - rna$Start + 1

# Special treatment for genes with negative start
rna$Length = ifelse(rna$Start > 0, rna$Length, rna$End - rna$Start)

# Calculate RPKM
rna$RPKM = 1000 * rna$RPM / rna$Length

# Combine RNAseq and ribosome profiling data
gen = rbind(
  rna[,grep("TotalReads", colnames(rna), invert=T)],
  rib[,grep("TotalReads", colnames(rib), invert=T)]
  )

# It is not interesting to include samples B, C, and D
gen = subset(gen, !(Sample %in% c("B","C","D")))

# Filter out low read count genes
# 3729 genes before
genes_below_threshold = unique(as.character(
  subset(gen, Reads < 32 & !(Sample %in% c("A","E")))$Name
  ))
gen = subset(gen, !(Name %in% genes_below_threshold))
# 2582 genes after

library(reshape2)

gen_wide = dcast(gen[,c("Sample","Name","RPKM")], Name ~ Sample, value.var="RPKM")

# Change column names
colnames(gen_wide)[grep("A", colnames(gen_wide), fixed=T)] = "RPF_RPKM_A"
colnames(gen_wide)[grep("E", colnames(gen_wide), fixed=T)] = "RPF_RPKM_E"
colnames(gen_wide)[grep("time0", colnames(gen_wide), fixed=T)] = "RNA_RPKM_A"
colnames(gen_wide)[grep("time24", colnames(gen_wide), fixed=T)] = "RNA_RPKM_E"

# Remove genes with no ribosome profiling data
gen_wide = subset(gen_wide, !is.na(RPF_RPKM_A) & !is.na(RPF_RPKM_E))

# TE
gen_wide$TE_A = gen_wide$RPF_RPKM_A / gen_wide$RNA_RPKM_A
gen_wide$TE_E = gen_wide$RPF_RPKM_E / gen_wide$RNA_RPKM_E

# log2 fold change ribosome profiling
gen_wide$log2_FC_RPF_RPKM = log2(gen_wide$RPF_RPKM_E / gen_wide$RPF_RPKM_A)

# log2 fold change RNAseq
gen_wide$log2_FC_RNA_RPKM = log2(gen_wide$RNA_RPKM_E / gen_wide$RNA_RPKM_A)

# log2 fold change translational efficiency
gen_wide$log2_FC_TE = log2(gen_wide$TE_E / gen_wide$TE_A)

# Add type and sequence
gen_wide = merge(gen_wide, gene_types)

# Reorder columns
# Name, Type, Sequence, RPF_RPKM_A, RPF_RPKM_E, RNA_RPKM_A, RNA_RPKM_E,
# log2_FC_RPF_RPKM, log2_FC_RNA_RPKM, TE_A, TE_E, log2_FC_TE

gen_wide = gen_wide[,c(
  "Name", "Type", "Sequence",
  "RPF_RPKM_A", "RPF_RPKM_E", "RNA_RPKM_A", "RNA_RPKM_E",
  "log2_FC_RPF_RPKM", "log2_FC_RNA_RPKM", "TE_A", "TE_E", "log2_FC_TE"
  )]

# Add annotations
gen_wide = merge(gen_wide, annotations, all.x = T)
# gen_wide$Process[is.na(gen_wide$Process)] = "No annotation"

# Save table
write.table(
  gen_wide, "analysis/TE_table.tab",
  sep = "\t", quote=F, row.names=F, col.names=T
)
