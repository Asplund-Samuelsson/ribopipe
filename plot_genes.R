#!/usr/bin/env Rscript


### COMMAND LINE ARGUMENTS #####################################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
indir = args[1] # A ribopipe results directory
genes = args[2] # Gene name(s); slr1834,sll1234, or interval; 8001-9001
outfile = args[3] # Output plot in PDF format

# TESTING
#indir="/ssd/jan/ribprof/seqdata/2017.03.29.PPE5/"
#genes="slr1329,sll1580,slr1908,sll1694"
#outfile="/tmp/ribopipe_plot.pdf"

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

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
rpm_strands = unlist(
  lapply(strsplit(rpm_files, "\\."), tail, n=1L)
)

# Obtain sample and strand IDs for readsPerGene files in order
gene_samples = unlist(
  lapply(lapply(strsplit(gene_files, "\\."), tail, n=3L), head, n=1L)
)
gene_strands = unlist(
  lapply(strsplit(gene_files, "\\."), tail, n=1L)
)


### LOAD DATA ##################################################################

# Load the files
rpm_data = lapply(rpm_files, read.delim, header=F)
gene_data = lapply(gene_files, read.delim, header=F)

# Add columns with sample and strand IDs to each dataframe
for( i in seq_along(rpm_data)){
  rpm_data[[i]] = cbind(rpm_data[[i]], Sample=rpm_samples[i])
  rpm_data[[i]] = cbind(rpm_data[[i]], Strand=rpm_strands[i])
}
for( i in seq_along(gene_data)){
  gene_data[[i]] = cbind(gene_data[[i]], Sample=gene_samples[i])
  gene_data[[i]] = cbind(gene_data[[i]], Strand=gene_strands[i])
}

# Create one dataframe for RPM data
rpm = do.call("rbind", rpm_data)
colnames(rpm)[1:2] = c("Position", "RPM")

# Create one dataframe for gene data
gen = do.call("rbind", gene_data)
colnames(gen)[1:4] = c("Name", "Start", "End", "Reads")


### PLOT DATA ##################################################################

suppressMessages(library(ggplot2))
suppressMessages(library(ggbio))
suppressMessages(library(GenomicRanges))

# Plotting for individual genes
gene_names = unlist(strsplit(genes, ","))

# Extract data for only the genes of interest
gdat = subset(gen, Name %in% gene_names)

# Expand gene data to every position from start to end
gdat_exp = data.frame(
  Name=character(),
  Position=numeric(),
  Sample=character(),
  Strand=character()
)

for( i in seq(1:(nrow(gdat)))){
  gdat_new = data.frame(
    Name=gdat[i,]$Name,
    Position=seq(gdat[i,]$Start, gdat[i,]$End),
    Sample=gdat[i,]$Sample,
    Strand=gdat[i,]$Strand
  )
  gdat_exp = rbind(gdat_exp, gdat_new)
}

# Merge with RPM data
grpm = merge(gdat_exp, rpm)

# Facilitate reversing order of positions in minus strand by multiplying with -1
grpm$Position = ifelse(grpm$Strand == "m", grpm$Position * -1, grpm$Position)

# Calculate minimum (=start) position per gene
mins = aggregate(Position ~ Name, grpm, min)
colnames(mins)[2] = "MinPosition"
grpm = merge(grpm, mins)

# Calculate gene positions
grpm$GenePosition = grpm$Position - grpm$MinPosition + 1

gp = ggplot(grpm, aes(x=GenePosition, y=RPM, fill=Strand))
gp = gp + geom_bar(position=position_dodge(), stat="identity")
gp = gp + theme_bw()
gp = gp + facet_grid(Sample~Name, scales="free_x", space="free_x")
gp = gp + theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
ggsave(outfile, gp, width=210/25.4, height=210/25.4)
