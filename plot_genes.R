#!/usr/bin/env Rscript


### COMMAND LINE ARGUMENTS #####################################################

# Load command line arguments
args = commandArgs(trailingOnly=T)
indir = args[1] # A ribopipe results directory
genes = args[2] # Gene name(s); slr1834,sll1234, or interval; ref_seq:8001:9001
plot_type = args[3] # Plot "facets" or "operon"
shift = as.numeric(args[4]) # Signal shift to apply to nucleotide RPM values
outfile = args[5] # Output plot in PDF format

# TESTING
# indir="/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick"
#
# genes="sll6098,slr1510,slr1511,sll1418"
# plot_type="facets"
# shift=-12
# outfile="/tmp/ribopipe_facets_plot.pdf"
#
# genes="slr1510,slr1511,sll1418"
# plot_type="operon"
# shift=-12
# outfile="/tmp/ribopipe_operon_plot.chromosome.pdf"
#
# genes="NC_005230.1:82000:87000"
# plot_type="operon"
# shift=-12
# outfile="/tmp/ribopipe_operon_plot.plasmid.pdf"
#
# genes="SGL_RS06570"
# plot_type="operon"
# shift=-12
# outfile="/tmp/ribopipe_operon_plot.SGL_RS06570.pdf"
#
# genes="SGL_RS14355"
# plot_type="operon"
# shift=-12
# outfile="/tmp/ribopipe_operon_plot.SGL_RS14355.pdf"
#
# genes="data/2018-05-15/plot_genes.test_batch.txt"
# plot_type="operon"
# shift=-12
# outfile="/tmp/ribopipe_operon_batch"
#
# genes="slr0611"
# shift=-12
# outfile="/tmp/slr0611.pdf"

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

message("Loading data...")

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

# Use data.table library for faster loading
suppressMessages(library(data.table))

# Load the files
rpm_data = lapply(rpm_files, fread, header=F)
gene_data = lapply(gene_files, fread, header=F)

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

### SHIFT RPM VALUES ###########################################################

library(dplyr)

if (shift) {

  message("Shifting RPM values...")

  # For the plus strand, the position shift is added
  # For the minus strand, the position shift is subtracted

  # Change position by the shift
  rpm_shift = rpm
  rpm_shift$Position = ifelse(
    rpm_shift$strand == "+",
    rpm_shift$Position + shift,
    rpm_shift$Position - shift
    )

  # Fold around beginning of circular genome
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

  rpm = rpm_shift[,c("Sequence","Position","RPM","Sample","strand")]

}

### DEFINE FUNCTIONS ###########################################################

# Load general libraries
suppressMessages(library(ggplot2))

# Expand gene data to every position
expand.genes = function(gdat){
  gdat_exp = data.frame(
    Sequence=character(),
    Name=character(),
    Position=numeric(),
    Sample=character(),
    strand=character()
    )

  for( i in seq(1:(nrow(gdat)))){
    gdat_new = data.frame(
      Sequence=gdat[i,]$Sequence,
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

### FACETS PLOTTING FUNCTION ###################################################

plot_facets = function(genes){

  message("Performing facets plotting...")

  # Plotting for individual genes
  gene_names = unlist(strsplit(genes, ","))

  # Extract data for only the genes of interest
  gdat = subset(gen, Name %in% gene_names)

  # Expand gene data to every position from start to end
  gdat_exp = expand.genes(gdat)

  # If any gene contains negative positions, give special treatment
  if(sum(gdat_exp$Position < 1)){

    if(!exists("sequence_sizes")){
      # Calculate size of sequences, if it wasn't done during shifting
      sequence_sizes = aggregate(Position ~ Sequence, rpm, max)
      colnames(sequence_sizes)[2] = "Size"
    }

    # Fold around beginning of circular genome
    gdat_exp$FakePosition = gdat_exp$Position
    gdat_exp = inner_join(gdat_exp, sequence_sizes)

    gdat_exp$Position = ifelse(
      gdat_exp$Position < 1,
      gdat_exp$Position + gdat_exp$Size,
      ifelse(
        gdat_exp$Position > gdat_exp$Size,
        gdat_exp$Position - gdat_exp$Size,
        gdat_exp$Position
        )
      )

    gdat_exp = gdat_exp[,
      c("Sequence","Name","Position","Sample","strand", "FakePosition")
    ]

  } else {
    gdat_exp$FakePosition = gdat_exp$Position
  }

  # Merge with RPM data
  grpm = inner_join(gdat_exp, filter(rpm, Position %in% unique(gdat_exp$Position)))

  # Facilitate reversing order of positions in minus strand by multiplying with -1
  grpm$Position = ifelse(grpm$strand == "-", grpm$FakePosition * -1, grpm$FakePosition)

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
  gp = gp + theme(
    axis.text.x = element_text(angle = 60, hjust=1, vjust=1),
    strip.background = element_blank()
  )
  gp = gp + colScale
  ggsave(outfile, gp, pdf, width=210/25.4, height=210/25.4)

}

### OPERON (GENOMIC RANGE) PLOTTING FUNCTION ###################################

plot_operon = function(genes){

  message("Performing operon plotting...")

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
    # Parse the specified range
    genes_list = unlist(strsplit(genes, ":"))
    sequence = genes_list[1]
    operon_start = as.numeric(genes_list[2])
    operon_end = as.numeric(genes_list[3])
    genes_full_range = c(operon_start, operon_end)
  }else{
    # Get range from specified genes
    gene_names = unlist(strsplit(genes, ","))
    gdat = subset(gen, Name %in% gene_names)
    genes_full_range = c(
      min(c(gdat$Start, gdat$End)) - 100,
      max(c(gdat$Start, gdat$End)) + 100
    )
    sequence = unique(gdat$Sequence)
    if (length(sequence) > 1) {
      stop("Attempting operon plot with multiple sequences.")
    }
  }

  # Calculate size of the sequence
  sequence_size = max(filter(rpm, Sequence == sequence)$Position)

  # Calculate wanted genome positions
  wanted_positions = (genes_full_range[1]):(genes_full_range[2])
  wanted_positions = ifelse(
    wanted_positions < 1,
    wanted_positions + sequence_size,
    ifelse(
      wanted_positions > sequence_size,
      wanted_positions - sequence_size,
      wanted_positions
    )
  )

  # Extract all gene data for the specified range
  gdat = subset(gen,
    (Start %in% wanted_positions | End %in% wanted_positions) &
    Sequence == sequence
  )

  # Translate end of sequence to negative values
  starts_at_end = gdat$Start > genes_full_range[2]
  gdat$Start = ifelse(
    starts_at_end,
    gdat$Start - sequence_size,
    gdat$Start
  )
  gdat$End = ifelse(
    starts_at_end,
    gdat$End - sequence_size,
    gdat$End
  )

  # Translate beginning of sequence to positive values
  ends_at_start = gdat$End < genes_full_range[1]
  gdat$Start = ifelse(
    ends_at_start,
    gdat$Start + sequence_size,
    gdat$Start
  )
  gdat$End = ifelse(
    ends_at_start,
    gdat$End + sequence_size,
    gdat$End
  )

  # Create gene names list if using range input
  if(grepl(":", genes, fixed=T)){gene_names = unique(gdat$Name)}

  # Clip genes to range
  gdat$Start[gdat$Start <= genes_full_range[1]] = genes_full_range[1]
  gdat$End[gdat$End >= genes_full_range[2]] = genes_full_range[2]

  # Expand gene data and merge with whole range of RPM data
  if(nrow(gdat)){
    gdat_exp = expand.genes(gdat)
    gdat_exp$Position = ifelse(
      gdat_exp$Position < 1,
      gdat_exp$Position + sequence_size,
      ifelse(
        gdat_exp$Position > sequence_size,
        gdat_exp$Position - sequence_size,
        gdat_exp$Position
        )
      )
    grpm = merge(
      gdat_exp,
      filter(rpm,
        Position %in% wanted_positions &
        Sequence == sequence
        ),
      all.y=T
      )
  }else{
    grpm = filter(rpm,
      Position %in% wanted_positions &
      Sequence == sequence
      )
  }

  # Add fake position
  grpm = inner_join(grpm, data.frame(
    Position = wanted_positions,
    FakePosition = (genes_full_range[1]):(genes_full_range[2])
  ))

  # Restore actual positions
  grpm$Position = grpm$FakePosition

  # Create GRanges object
  gene_ranges = unique(gdat[,c("Name", "Start", "End", "strand")])

  # Check the dominant direction of the genes
  x_reverse = names(sort(
      table(subset(gene_ranges, Name %in% gene_names)$strand),
      decreasing=TRUE
  )[1]) == "-"

  # If there is no dominant direction, set x_reverse to FALSE
  no_dom_dir = length(unique(as.data.frame(
    table(subset(gene_ranges, Name %in% gene_names)$strand)
  )$Freq)) == 1
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
    if (1 %in% plot_data$Position){
      gp = gp + geom_vline(
        xintercept=1, alpha=0.8,
        colour="red", linetype="dashed", size=0.2
      )
    }
    if ((sequence_size + 1) %in% plot_data$Position){
      gp = gp + geom_vline(
        xintercept=(sequence_size + 1), alpha=0.8,
        colour="red", linetype="dashed", size=0.2
      )
    }
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
          xintercept=position, colour="black", size=0.2, alpha=0.4
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
    gp = gp + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank()
    )
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
  garbage = dev.off()

}

### PERFORM PLOTTING ###########################################################

if (file.exists(genes)){
  # If the "genes" variable is a batch file, plot data for each line
  lines = scan(genes, character())
} else {
  # If the "genes" variable is a single line, plot data for only that
  lines = genes
}

original_outfile = outfile

# Perform the plotting
for (genes in lines){
  # Outfile is not complete if in batch mode
  if (length(lines) > 1){
    outfile = paste(
      c(original_outfile, ".", gsub(",|:", "_", genes), ".pdf"), collapse=""
    )
  }
  # Plot facets
  if (tolower(plot_type) == "facets"){
    plot_facets(genes)
  }
  # Plot operon
  if (tolower(plot_type) == "operon"){
    plot_operon(genes)
  }
}

message("Done.")
