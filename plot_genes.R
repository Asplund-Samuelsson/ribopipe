#!/usr/bin/env Rscript


### COMMAND LINE ARGUMENTS #####################################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
indir = args[1] # A ribopipe results directory
genes = args[2] # Gene name(s); slr1834,sll1234, or interval; 8001:9001
plot_type = args[3] # Plot "facets" or "operon"
outfile = args[4] # Output plot in PDF format

# TESTING
#indir="/ssd/jan/ribprof/seqdata/2017.03.29.PPE5/"
#genes="sll1577,sll1578,sll1579,sll1580,ssl3093"
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


### LOAD LIBRARIES #############################################################

suppressMessages(library(ggplot2))
suppressMessages(library(ggbio))
suppressMessages(library(GenomicRanges))
library(egg)


### DEFINE FUNCTIONS ###########################################################

# Expand gene data to every position
expand.genes = function(gdat){
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

  return(gdat_exp)
}


### PLOT DATA AS FACETS ########################################################

if (tolower(plot_type) == "facets"){

  # Plotting for individual genes
  gene_names = unlist(strsplit(genes, ","))

  # Extract data for only the genes of interest
  gdat = subset(gen, Name %in% gene_names)

  # Expand gene data to every position from start to end
  gdat_exp = expand.genes(gdat)

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

  # Make Name factor for preserving user-specified order
  grpm$Name = factor(as.character(grpm$Name), levels=gene_names)

  gp = ggplot(grpm, aes(x=GenePosition, y=RPM, fill=Strand))
  gp = gp + geom_bar(position=position_dodge(), stat="identity")
  gp = gp + theme_bw()
  gp = gp + facet_grid(Sample~Name, scales="free_x", space="free_x")
  gp = gp + theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
  ggsave(outfile, gp, width=210/25.4, height=210/25.4)

}

### PLOT GENOMIC RANGE (OPERON) ################################################

if (tolower(plot_type) == "operon"){

  # Check if range or set of gene names
  if(grepl(":", genes, fixed=T)){
    # Get genes from specified range
    genes_full_range = as.numeric(unlist(strsplit(genes, ":")))
    gdat = subset(gen, Start >= genes_full_range[1] & End <= genes_full_range[2])
  }else{
    # Get range from specified genes
    gene_names = unlist(strsplit(genes, ","))
    gdat = subset(gen, Name %in% gene_names)
    genes_full_range = c(
      min(c(gdat$Start, gdat$End)) - 100,
      max(c(gdat$Start, gdat$End)) + 100
    )
  }

  # Expand gene data and merge with whole range of RPM data
  gdat_exp = expand.genes(gdat)
  grpm = merge(
    gdat_exp,
    subset(rpm, Position >= genes_full_range[1] & Position <= genes_full_range[2]),
    all.y=T
    )

  # Create GRanges object
  gene_ranges = unique(gdat[,c("Name", "Start", "End", "Strand")])
  gene_ranges$Middle = (gene_ranges$Start + gene_ranges$End)/2
  gene_ranges$Y = 1
  gene_ranges$LabelPosition = ifelse(
    gene_ranges$Strand == "p",
    gene_ranges$Start + 10,
    gene_ranges$End - 10
    )
  gene_ranges$hjust = ifelse(gene_ranges$Strand == "p", "right", "left")
  gr = GRanges(
    seqnames = genes,
    ranges = IRanges(
      start = gene_ranges$Start,
      width = gene_ranges$End - gene_ranges$Start + 1
    ),
    strand = ifelse(gene_ranges$Strand == "p", "+", "-"),
    Name = gene_ranges$Name
    )

  # Perform plotting

  # Plot RPM
  gp = ggplot(grpm, aes(x=Position, y=RPM))
  gp = gp + geom_bar(position=position_dodge(), stat="identity")
  gp = gp + theme_bw()
  gp = gp + facet_grid(Sample~.)
  if(names(sort(table(gene_ranges$Strand),decreasing=TRUE)[1]) == "m"){
    gp = gp + scale_x_reverse()
  }
  gp = gp + theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())
  gp = gp + xlim(genes_full_range)
  gp_rpm = gp

  # Plot genes
  gp = ggplot(gr)
  gp = gp + geom_hline(yintercept=1)
  gp = gp + geom_alignment(
    range.geom = "arrowrect",
    gap.geom = "segment",
    length = 1
    )
  if(names(sort(table(gene_ranges$Strand),decreasing=TRUE)[1]) == "m"){
    gp = gp + scale_x_reverse()
  }
  gp = gp + geom_text(
    aes(x=LabelPosition, y=Y, label=Name, hjust=hjust), gene_ranges,
    colour="white", size=2.5
    )
  gp = gp + theme_bw()
  gp = gp + xlim(genes_full_range)
  gp_genes = gp@ggplot # @ggplot FOR PLOTTING WITH REGULAR GGPLOT OBJECTS

  # Align plots
  pdf(outfile, width=210/25.4, height=140/25.4, onefile=FALSE)
  ggarrange(gp_rpm, gp_genes, ncol=1, heights=c(20,1))
  dev.off()

}
