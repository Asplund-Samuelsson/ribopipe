
### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input directories
indir="/ssd/common/proj/RibosomeProfiling/ribopipe_testing/CSD2_3prime"

# Define codon file
codon_file = "data/2017-11-20/Synechocystis_PCC_6803_chr.codon_seqpos_strand_genpos123.tab"

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
codons = read.table(codon_file, header=F, sep="\t", stringsAsFactors=F)

# Update codon header
colnames(codons) = c("Codon", "SeqPos", "Strand", "GenPos1", "GenPos2", "GenPos3")

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

### NORMALIZE ##################################################################

# Count total reads mapped to genes
total_count = aggregate(Reads ~ Sample, gen, sum)
colnames(total_count)[2] = "TotalReads"

# Merge with gen to create new column with total read count per sample
gen = merge(gen, total_count)

# Calculate RPM
gen$RPM = gen$Reads * 1000000 / gen$TotalReads

### ANALYZE ####################################################################

# Give codons a unique identifier
codons$CodonID = rownames(codons)

library(reshape2)
codons_wide = melt(codons, id.vars=c("CodonID", "Codon", "SeqPos", "Strand"))
colnames(codons_wide)[c(4,6)] = c("strand", "Position")

# Add gene names

# Use the first nucleotide of the codon as its position for gene matching
codon_positions = unique(subset(codons_wide, variable == "GenPos1")[,c("CodonID","Position","strand")])

# Replace negative position in genome with proper position
gen$Start[gen$Start == -200] = max(codons_wide$Position) - 200 + 1

# Match codons with the corresponding genes via position
gene_start_end = unique(gen[,c("Name","strand","Start","End")])

library(parallel)
cl = makeCluster(30)
clusterExport(cl, "gene_start_end")

codon_gene_name = parApply(
  cl, codon_positions, 1,
  function (x) {
    Pos = as.numeric(x[2])
    Strand = as.character(x[3])
    S = subset(gene_start_end, Start <= Pos & End >= Pos & strand == Strand)
    as.character(S$Name)
  }
)

stopCluster(cl)

# Replace the empty character variables (character(0)) with an empty string
codon_gene_name[unlist(lapply(codon_gene_name, function(x){identical(x, character(0))}))] = ""

# Replace multiple gene position name vectors in list with "MULTI"
codon_gene_name[unlist(lapply(codon_gene_name, length)) > 1] = "MULTI"

# Create the vector...
codon_gene_name_vec = unlist(codon_gene_name)
# Length is finally correct

# Add Name column
codon_positions$Name = codon_gene_name_vec

# Fix the odd gene that straddles the start of the chromosome
codon_positions$Name[codon_positions$CodonID %in% 1:324] = "slr0611"

# Filter out MULTI codons
codon_positions = subset(codon_positions, Name != "MULTI")

# Calculate gene length
gen$Length = gen$End - gen$Start + 1

# Special treatment for slr0611
gen$Length[gen$Name == "slr0611"] = 772 + 200

# Calculate RPKM
gen$RPKM = 1000 * gen$RPM / gen$Length

# Extend the RPM table with shifted values from position -10 to +40
# For the plus strand, the position shift is added
# For the minus strand, the position shift is subtracted

# Merge codons with RPM values for all positions and shifts

cl = makeCluster(10) # USE NO MORE THAN 10 CLUSTERS WITH 128 GB RAM
clusterExport(cl, c("rpm", "codons_wide", "gen", "codon_positions"))

codon_rpm_list = parLapply(
    cl, -50:50, function(x){

      # Change position by the shift
      rpm_shift = rpm
      rpm_shift$Shift = x
      rpm_shift$Position = ifelse(
        rpm_shift$strand == "+",
        rpm_shift$Position + rpm_shift$Shift,
        rpm_shift$Position - rpm_shift$Shift
        )

      # Fold around beginning of circular genome
      rpm_shift$Position = ifelse(
        rpm_shift$Position < 1,
        rpm_shift$Position + max(rpm$Position),
        ifelse(
          rpm_shift$Position > max(rpm$Position),
          rpm_shift$Position - max(rpm$Position),
          rpm_shift$Position
          )
        )

      # Calculate the mean RPM shift per codon
      codon_rpm = aggregate(
        RPM ~ strand + CodonID + Codon + Sample + Shift,
        merge(codons_wide, rpm_shift),
        mean
        )

      # Add gene name
      codon_rpm = merge(codon_rpm, codon_positions[,c("CodonID","Name")])

      # Add gene RPKM to codon_rpm
      codon_rpm = merge(codon_rpm, gen[,c("Sample","strand","Name","RPKM")])

      # Normalize RPM codon to RPKM of gene
      codon_rpm$PauseScore = 1000 * codon_rpm$RPM / codon_rpm$RPKM

      # Remove instances where PauseScore is NaN or infinite
      codon_rpm = codon_rpm[is.finite(codon_rpm$PauseScore),]

      # Calculate median PauseScore for all codons and samples
      codon_median = aggregate(PauseScore ~ Codon + Sample, codon_rpm, median)
      colnames(codon_median)[3] = "MedianPauseScore"

      # Calculate PauseScore MAD for all codons and samples
      codon_MAD = aggregate(PauseScore ~ Codon + Sample, codon_rpm, mad)
      colnames(codon_MAD)[3] = "PauseScoreMAD"

      # Return results
      codon_rpm = merge(codon_rpm, merge(codon_median, codon_MAD))
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
  "results/2018-01-12/codon_PauseScore_3prime_with_shift.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )

### PLOT DATA ##################################################################

# Create a table with the median codon pause scores per sample and shift
psMAD = unique(codon_rpm[,c("Codon", "Sample", "Shift", "MedianPauseScore", "PauseScoreMAD")])
colnames(psMAD) = c("Codon", "Sample", "Shift", "PauseScore", "psMAD")

# Save the table
write.table(
  psMAD,
  "results/2018-01-12/codon_PauseScore_medians_3prime_with_shift.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )

# Add colours for amino acids
colour_file = "../main/art/2016-09-30/colours.txt"
colours = scan(colour_file, character(), quote = "")
colours = c("#bdbdbd", colours, "#fc4e2a") # Add additional AA and STOP colours

# Load codon translation table
translation_file = "data/2017-11-21/translation_table_11.tab"
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
gp = gp + geom_vline(xintercept=-14, color="black")
gp = gp + geom_errorbar()
gp = gp + geom_line()
gp = gp + geom_point()
gp = gp + theme_bw()
gp = gp + scale_x_continuous(breaks=-50:50)
gp = gp + facet_grid(Sample ~ Codon + AminoAcid)
gp = gp + scale_colour_manual(values=colours)
gp = gp + theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=7))

ggsave("art/2018-01-12/pause_score_distribution.codon_median.3prime_shift.1.pdf", gp, width=14000/25.4, height=180/25.4, limitsize=F)
