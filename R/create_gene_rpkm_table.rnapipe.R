#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
RNAdir = "."
genelist_file = args[1]

# # Load command arguments manually
# RNAdir = "/hdd/common/proj/RNAseq/results/2018-09-14_CCE1/2018.06.11.CCE1_RNAseq"
# genelist_file = "/ssd/common/tools/ribopipe/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"

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

# Calculate the total number of reads on all genes in each sample
tc_rna = merge(tc_rna, aggregate(Reads~Sample, rna, sum))
colnames(tc_rna)[3] = "totalReads_genes"

# Calculate the total number of reads on CDSs in each sample
tc_rna = merge(
  tc_rna,
  aggregate(Reads~Sample,
            subset(rna, Name %in% subset(gene_types, Type == "CDS")$Name),
            sum))
colnames(tc_rna)[4] = "totalReads_orf"

# Add total count to per-gene RPM data
rna = merge(rna, tc_rna)

# Calculate RPM
rna$RPM = rna$Reads * 1000000 / rna$totalReads_genes

# Calculate gene length
rna$Length = rna$End - rna$Start + 1

# Special treatment for genes with negative start
rna$Length = ifelse(rna$Start > 0, rna$Length, rna$End - rna$Start)

# Calculate RPKM
rna$RPKM = 1000 * rna$RPM / rna$Length

# Subset to coding sequences only
rna_CDS = subset(rna, Name %in% subset(gene_types, Type == "CDS")$Name)

# Recalculate RPM and RPKM. scalingfactor = totalReads_orf. 
rna_CDS$RPM = rna_CDS$RPM * rna_CDS$totalReads_genes /
  rna_CDS$totalReads_orf
rna_CDS$RPKM = rna_CDS$RPKM * rna_CDS$totalReads_genes /
  rna_CDS$totalReads_orf

### WRITE TABLE ################################################################

write.table(
  rna_CDS,
  "analysis/all_CDS_RPKM_no_filter.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )

write.table(
  rna,
  "analysis/all_genes_RPKM_no_filter.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )
