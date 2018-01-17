#!/usr/bin/env Rscript


### COMMAND LINE ARGUMENTS #####################################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
indir = args[1] # A ribopipe results directory
genes = args[2] # Gene name(s); slr1834,sll1234, or interval; 8001:9001
plot_type = args[3] # Plot "facets" or "operon"
outfile = args[4] # Output plot in PDF format

# TESTING
#indir="/ssd/common/proj/RibosomeProfiling/ribopipe_testing/CSD2_3prime"
#genes="slr1510,slr1511,sll1418"
#outfile="/tmp/ribopipe_plot.pdf"
#plot_type="operon"

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

### LOAD DATA ##################################################################

# Load the files
rpm_data = lapply(rpm_files, read.delim, header=F)
gene_data = lapply(gene_files, read.delim, header=F)

# Add columns with sample and strand IDs to each dataframe
for( i in seq_along(rpm_data)){
  rpm_data[[i]] = cbind(rpm_data[[i]], Sample=rpm_samples[i])
  rpm_data[[i]] = cbind(rpm_data[[i]], strand=rpm_strands[i])
}
for( i in seq_along(gene_data)){
  gene_data[[i]] = cbind(gene_data[[i]], Sample=gene_samples[i])
  gene_data[[i]] = cbind(gene_data[[i]], strand=gene_strands[i])
}

# Create one dataframe for RPM data
rpm = do.call("rbind", rpm_data)
colnames(rpm)[1:2] = c("Position", "RPM")

# Create one dataframe for gene data
gen = do.call("rbind", gene_data)
colnames(gen)[1:4] = c("Name", "Start", "End", "Reads")

### DEFINE FUNCTIONS ###########################################################

# Load general libraries
library(ggplot2)

# Expand gene data to every position
expand.genes = function(gdat){
  gdat_exp = data.frame(
    Name=character(),
    Position=numeric(),
    Sample=character(),
    strand=character()
    )

  for( i in seq(1:(nrow(gdat)))){
    gdat_new = data.frame(
      Name=gdat[i,]$Name,
      Position=seq(gdat[i,]$Start, gdat[i,]$End),
      Sample=gdat[i,]$Sample,
      strand=gdat[i,]$strand
      )
    gdat_exp = rbind(gdat_exp, gdat_new)
  }

  return(gdat_exp)
}

# Define strand color scale
strandColors = c("#b35806", "#542788")
names(strandColors) = c("+","-")
colScale <- scale_fill_manual(name = "strand", values = strandColors)

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
  grpm$Position = ifelse(grpm$strand == "-", grpm$Position * -1, grpm$Position)

  # Calculate minimum (=start) position per gene
  mins = aggregate(Position ~ Name, grpm, min)
  colnames(mins)[2] = "MinPosition"
  grpm = merge(grpm, mins)

  # Calculate gene positions
  grpm$GenePosition = grpm$Position - grpm$MinPosition + 1

  # Make Name factor for preserving user-specified order
  grpm$Name = factor(as.character(grpm$Name), levels=gene_names)

  gp = ggplot(grpm, aes(x=GenePosition, y=RPM, fill=strand))
  gp = gp + geom_bar(position=position_dodge(), stat="identity")
  gp = gp + theme_bw()
  gp = gp + facet_grid(Sample~Name, scales="free_x", space="free_x")
  gp = gp + scale_x_continuous(
    breaks = seq(0, max(grpm$GenePosition), 200),
    minor_breaks = seq(0, max(grpm$GenePosition), 20)
    )
  gp = gp + theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))
  gp = gp + colScale
  ggsave(outfile, gp, pdf, width=210/25.4, height=210/25.4)

}

### PLOT GENOMIC RANGE (OPERON) ################################################

if (tolower(plot_type) == "operon"){

  # Load libraries
  suppressMessages(library(ggbio))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(egg))
  suppressMessages(library(ggrepel))
  suppressMessages(library(plyr))

  # Clean up when testing
  if(exists("gr")){rm("gr")}

  # Check if range or set of gene names
  if(grepl(":", genes, fixed=T)){
    # Get genes from specified range
    genes_full_range = as.numeric(unlist(strsplit(genes, ":")))
    gdat = subset(gen, (Start >= genes_full_range[1] & Start <= genes_full_range[2]) | (End >= genes_full_range[1] & End <= genes_full_range[2]))
    # Clip genes to range
    gdat$Start[gdat$Start <= genes_full_range[1]] = genes_full_range[1]
    gdat$End[gdat$End >= genes_full_range[2]] = genes_full_range[2]
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
  if(nrow(gdat)){
    gdat_exp = expand.genes(gdat)
    grpm = merge(
      gdat_exp,
      subset(rpm, Position >= genes_full_range[1] & Position <= genes_full_range[2]),
      all.y=T
      )
  }else{
    grpm = subset(rpm, Position >= genes_full_range[1] & Position <= genes_full_range[2])
  }

  # Create GRanges object
  gene_ranges = unique(gdat[,c("Name", "Start", "End", "strand")])

  # Check the dominant direction of the genes
  x_reverse = names(sort(table(gene_ranges$strand),decreasing=TRUE)[1]) == "-"

  # If there is no dominant direction, set x_reverse to FALSE
  no_dom_dir = length(unique(as.data.frame(table(gene_ranges$strand))$Freq)) == 1
  if(no_dom_dir){
    x_reverse = FALSE
  }

  # Add data for plot adjustment
  if (nrow(gene_ranges)){
    gene_ranges$Middle = (gene_ranges$Start + gene_ranges$End)/2
    gene_ranges$Y = 1
    gene_ranges$LabelPosition = ifelse(
      gene_ranges$strand == "+",
      gene_ranges$Start + 10,
      gene_ranges$End - 10
      )
    if (x_reverse){
      gene_ranges$hjust = ifelse(gene_ranges$strand == "+", "right", "left")
    }else{
      gene_ranges$hjust = ifelse(gene_ranges$strand == "+", "left", "right")
    }
    gr = GRanges(
      seqnames = genes,
      ranges = IRanges(
        start = gene_ranges$Start,
        width = gene_ranges$End - gene_ranges$Start + 1
      ),
      strand = gene_ranges$strand,
      Name = gene_ranges$Name
      )
  }

  # Perform plotting

  # Define breaks
  x_major_breaks_magnitude = 10^round_any(
    log10(genes_full_range[2] - genes_full_range[1] + 1), 1, f=floor
    )
  x_minor_breaks_magnitude = x_major_breaks_magnitude / 10

  x_major_breaks = seq(
    round_any(genes_full_range[1], x_major_breaks_magnitude, f=ceiling),
    round_any(genes_full_range[2], x_major_breaks_magnitude, f=floor),
    x_major_breaks_magnitude
    )
  x_minor_breaks = seq(
    round_any(genes_full_range[1], x_minor_breaks_magnitude, f=ceiling),
    round_any(genes_full_range[2], x_minor_breaks_magnitude, f=floor),
    x_minor_breaks_magnitude
    )

  # Plot RPM
  plot_RPM = function(S, reverse_facets){
    # Subset to selected strand
    plot_data = subset(grpm, strand == S)
    # Reorder samples if the facet order should be reversed
    samples = unique(as.character(plot_data$Sample))
    samples = samples[order(samples)]
    if (reverse_facets){
      samples = samples[order(samples, decreasing=T)]
    }
    plot_data$Sample = factor(as.character(plot_data$Sample), levels=samples)
    gp = ggplot(plot_data, aes(x=Position, y=RPM, fill=strand))
    if (S == "+"){
      gp = gp + geom_bar(position=position_dodge(), stat="identity")
    }else{
      gp = gp + geom_bar(position=position_dodge(), stat="identity")
    }
    if(exists("gr")){
      # If there is a GenomicRanges object, add vlines indicating gene bounds
      gr_sub = subset(gr, strand == S)
      gr_df = as.data.frame(gr_sub@ranges)
      gene_bounds = c(gr_df$start, gr_df$end)
      for (position in gene_bounds){
        gp = gp + geom_vline(
          xintercept=position, colour="black", size=0.2, alpha=0.8
          )
      }
    }
    gp = gp + theme_bw()
    gp = gp + facet_grid(Sample~.)
    # Set up x axis
    if(x_reverse){
      gp = gp + scale_x_reverse(
        limits=rev(genes_full_range), expand=c(0,0),
        breaks=x_major_breaks, minor_breaks=x_minor_breaks
        )
    }
    if(!x_reverse){
      gp = gp + scale_x_continuous(
        limits=genes_full_range, expand=c(0,0),
        breaks=x_major_breaks, minor_breaks=x_minor_breaks
        )
    }
    # Set up y axis
    y_axis_limits = c(0, max(grpm$RPM)*1.05)
    if((!x_reverse & S == "-") | (x_reverse & S == "+")){
      gp = gp + scale_y_reverse(limits=rev(y_axis_limits))
    }else{
      gp = gp + scale_y_continuous(limits=y_axis_limits)
    }
    gp = gp + colScale
    # Remove x axis decorations
    gp = gp + theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank())
    return(gp)
  }

  # Plot genes
  if(exists("gr")){
    gp = ggplot(gr, aes(fill=strand))
    gp = gp + geom_hline(yintercept=1)
    gp = gp + geom_alignment(
      range.geom = "rect",
      gap.geom = "segment"
      )
    if(x_reverse){
      gp = gp + scale_x_reverse(
        expand=c(0,0),
        breaks=x_major_breaks, minor_breaks=x_minor_breaks
        )
    }else{
      gp = gp + scale_x_continuous(
        expand=c(0,0),
        breaks=x_major_breaks, minor_breaks=x_minor_breaks
        )
    }
    gp = gp + geom_label_repel(
      aes(x=LabelPosition, y=Y, label=Name), gene_ranges,
      colour="lightgrey", size=2.5, min.segment.length=unit(0, "lines"),
      alpha=0.8
      )
    gp = gp + theme_bw()
    if(x_reverse){
      gp = gp + xlim(rev(genes_full_range))
    }else{
      gp = gp + xlim(genes_full_range)
    }
    gp = gp + scale_y_continuous(expand=c(0.2,0.2))
    gp = gp + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
    gp = gp + colScale
    gp = gp + guides(fill=F, colour=F)
    gp_genes = gp@ggplot # @ggplot FOR PLOTTING WITH REGULAR GGPLOT OBJECTS
  }else{
    dummy_gr_data = data.frame(
      Position=seq(genes_full_range[1], genes_full_range[2]),
      Y=1
      )
    gp = ggplot(dummy_gr_data, aes(x=Position, y=Y))
    gp = gp + geom_line()
    gp = gp + geom_hline(yintercept=1)
    gp = gp + theme_bw()
    if(x_reverse){
      gp = gp + xlim(rev(genes_full_range))
    }else{
      gp = gp + xlim(genes_full_range)
    }
    gp = gp + theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
    gp_genes = gp
  }

  # Align plots
  pdf(outfile, width=210/25.4, height=270/25.4, onefile=FALSE)
  if(!no_dom_dir & x_reverse){
    gp_rpm_p = plot_RPM("+", reverse_facets=T)
    gp_rpm_m = plot_RPM("-", reverse_facets=F)
    ggarrange(gp_rpm_m, gp_genes, gp_rpm_p, ncol=1, heights=c(20,2,20))
  }else{
    gp_rpm_p = plot_RPM("+", reverse_facets=F)
    gp_rpm_m = plot_RPM("-", reverse_facets=T)
    ggarrange(gp_rpm_p, gp_genes, gp_rpm_m, ncol=1, heights=c(20,2,20))
  }
  dev.off()

}
