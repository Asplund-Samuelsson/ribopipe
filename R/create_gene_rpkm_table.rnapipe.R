#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
RNAdir = "."
genelist_file = args[1]

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

# Load gene types
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

# Subset to coding sequences only
rna = subset(rna, Name %in% subset(gene_types, Type == "CDS")$Name)

### WRITE TABLE ################################################################

write.table(
  rna,
  "analysis/all_CDS_RPKM_no_filter.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )
