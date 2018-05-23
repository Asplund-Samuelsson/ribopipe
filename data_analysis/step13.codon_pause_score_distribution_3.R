#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# translation_file = "tools/ribopipe/data_analysis/translation_table_11.tab"
# colour_file = "tools/ribopipe/data_analysis/colours.txt"
# indir = "/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick"
# codon_file = "data/2018-03-28/Syn_PCC6803.1_chro_7plasmids.seq_gene_codon_seqpos_strand_genpos123.tab"

# Load command line arguments
args = commandArgs(trailingOnly=T)
translation_file = args[1]
colour_file = args[2]

# Define input directories
indir="."

# Define codon file
codon_file = "analysis/codon_seqpos.tab"

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
library(dplyr)

# Load the files
rpm_data = lapply(rpm_files, fread, header=F)
gene_data = lapply(gene_files, fread, header=F)
codons = read.table(codon_file, header=F, sep="\t", stringsAsFactors=F)
gen_tc_data = unlist(lapply(gen_tc_files, scan))

# Update codon header
colnames(codons) = c(
  "Sequence", "Name", "Codon", "SeqPos",
  "Strand", "GenPos1", "GenPos2", "GenPos3"
  )

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
gen = subset(gen, Name %in% codons$Name)

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

# Remove all genes that have less than 128 reads in any sample
genes_below_threshold = as.character(unique(subset(gen, Reads < 128)$Name))

gen = subset(gen, !(Name %in% genes_below_threshold))

# Calculate gene length
gen$Length = gen$End - gen$Start + 1

# Special treatment for genes with negative start
gen$Length = ifelse(gen$Start > 0, gen$Length, gen$End - gen$Start)

# Calculate RPKM
gen$RPKM = 1000 * gen$RPM / gen$Length

### ANALYZE ####################################################################

# Give codons a unique identifier
codons$CodonID = rownames(codons)

library(reshape2)
codons_wide = melt(
  codons,
  id.vars=c("Sequence", "Name", "CodonID", "Codon", "SeqPos", "Strand")
  )
colnames(codons_wide)[c(6,8)] = c("strand", "Position")

# Extend the RPM table with shifted values from position -10 to +40
# For the plus strand, the position shift is added
# For the minus strand, the position shift is subtracted

# Merge codons with RPM values for all positions and shifts

library(parallel)

cl = makeCluster(5) # USE NO MORE THAN 10 CLUSTERS WITH 128 GB RAM
clusterExport(
  cl,
  c("rpm", "codons_wide", "gen", "sequence_sizes", "tc_gen")
  )

codon_rpm_list = parLapply(
    cl, -50:50, function(x){

      library(dplyr)

      # Change position by the shift
      rpm_shift = rpm
      rpm_shift$Shift = x
      rpm_shift$Position = ifelse(
        rpm_shift$strand == "+",
        rpm_shift$Position + rpm_shift$Shift,
        rpm_shift$Position - rpm_shift$Shift
        )

      # Fold around beginning of circular genome
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

      # Recalculate read count
      rpm_shift = inner_join(rpm_shift, tc_gen)
      rpm_shift$ReadsNuc = rpm_shift$RPM * rpm_shift$TotalReads / 1000000

      # Calculate the sum of shifted reads per codon
      codon_rpm = aggregate(
        ReadsNuc ~ Sequence + strand + CodonID + Codon + Sample + Shift,
        inner_join(codons_wide, rpm_shift),
        sum
        )

      # Add gene name
      codon_rpm = inner_join(codon_rpm, unique(codons_wide[,c("CodonID","Name")]))

      # Add gene Reads to codon_rpm
      codon_rpm = inner_join(
        codon_rpm,
        gen[,c("Sample", "strand", "Name", "Reads", "Length")]
        )

      # Calculate Pause Score
      codon_rpm$PauseScore = (codon_rpm$ReadsNuc / 3) / (codon_rpm$Reads / codon_rpm$Length)

      # Remove instances where PauseScore is NaN or infinite
      codon_rpm = codon_rpm[is.finite(codon_rpm$PauseScore),]

      # Calculate median PauseScore for all codons and samples
      codon_median = aggregate(PauseScore ~ Codon + Sample, codon_rpm, median)
      colnames(codon_median)[3] = "MedianPauseScore"

      # Calculate PauseScore MAD for all codons and samples
      codon_MAD = aggregate(PauseScore ~ Codon + Sample, codon_rpm, mad)
      colnames(codon_MAD)[3] = "PauseScoreMAD"

      # Return results
      codon_rpm = inner_join(codon_rpm, inner_join(codon_median, codon_MAD))
      codon_rpm

    }
  )

stopCluster(cl)

# Combine into one data frame
library(data.table)
codon_rpm = as.data.frame(rbindlist(codon_rpm_list))

# Save a table
write.table(
  codon_rpm,
  "analysis/codon_PauseScore.all_shifts.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )

### PLOT DATA ##################################################################

# Create a table with the median codon pause scores per sample and shift
psMAD = unique(codon_rpm[,c("Codon", "Sample", "Shift", "MedianPauseScore", "PauseScoreMAD")])
colnames(psMAD) = c("Codon", "Sample", "Shift", "PauseScore", "psMAD")

# Save the table
write.table(
  psMAD,
  "analysis/codon_PauseScore_medians.all_shifts.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )

# Add colours for amino acids
colours = scan(colour_file, character(), quote = "")
colours = c("#bdbdbd", colours, "#fc4e2a") # Add additional AA and STOP colours

# Load codon translation table
t_table = read.table(translation_file, header=F, sep="\t", stringsAsFactors=F)
colnames(t_table) = c("Codon", "AminoAcid")
amino_acids = c(grep("Ter", unique(t_table$AminoAcid), value=T, invert=T), "Ter")

# Add amino acids to pause score data
psMAD = merge(psMAD, t_table)
psMAD$AminoAcid = factor(psMAD$AminoAcid, levels=amino_acids)

# Sort Codons by median pausescore
psMAD_median_means = aggregate(PauseScore ~ Codon, psMAD, mean)
psMAD_median_means = psMAD_median_means[order(psMAD_median_means$PauseScore),]
psMAD$Codon = factor(psMAD$Codon, levels=psMAD_median_means$Codon)

# Plot it
library(ggplot2)

# Calculate max and min values
psMAD$Max = psMAD$PauseScore + psMAD$psMAD
psMAD$Min = psMAD$PauseScore - psMAD$psMAD

psMAD = merge(psMAD, t_table)
psMAD$AminoAcid = factor(psMAD$AminoAcid, levels=amino_acids)

# Sort Codons by median pausescore
psMAD_median_means = aggregate(PauseScore ~ Codon, psMAD, mean)
psMAD_median_means = psMAD_median_means[order(psMAD_median_means$PauseScore),]
psMAD$Codon = factor(psMAD$Codon, levels=psMAD_median_means$Codon)

gp = ggplot(psMAD, aes(x=Shift, y=PauseScore, ymax=Max, ymin=Min, colour=AminoAcid))
gp = gp + geom_vline(xintercept=-12, color="black")
gp = gp + geom_errorbar()
gp = gp + geom_line()
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + scale_x_continuous(breaks=-50:50)
gp = gp + facet_grid(Sample ~ Codon + AminoAcid)
gp = gp + scale_colour_manual(values=colours)
gp = gp + theme(
  axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7),
  strip.background = element_blank()
)

ggsave("analysis/pause_score_distribution.codon_median.all_shifts.pdf", gp, width=14000/25.4, height=180/25.4, limitsize=F)
