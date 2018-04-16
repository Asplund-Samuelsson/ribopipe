
### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input directories
indir="/hdd/common/proj/RibosomeProfiling/results/2018-03-26/CSD2_plasmid_test_2"

# Define codon file
codon_file = "data/2018-03-28/Syn_PCC6803.1_chro_7plasmids.seq_gene_codon_seqpos_strand_genpos123.tab"

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
  codons, id.vars=c("CodonID", "Codon", "Sequence", "Name", "SeqPos", "Strand")
  )
colnames(codons_wide)[c(6,8)] = c("strand", "Position")

# Calculate the mean RPM shift per codon
codon_rpm = aggregate(
  RPM ~ strand + CodonID + Codon + Sample + Sequence + Name,
  inner_join(
    codons_wide,
    gen_allpos[,c("Sequence", "Name", "Sample", "strand", "Position", "RPM")]
    ),
  mean
  )

# Add gene RPKM to codon_rpm
codon_rpm = inner_join(codon_rpm, gen[,c("Sample","Name","RPKM")])

# Normalize RPM codon to RPKM of gene
codon_rpm$PauseScore = 1000 * codon_rpm$RPM / codon_rpm$RPKM

# Remove instances where PauseScore is NaN or infinite
codon_rpm = codon_rpm[is.finite(codon_rpm$PauseScore),]


### CALCULATE BOXPLOT VALUES ###################################################

# Calculate PauseScore Quantiles for all codons and samples
codon_quantiles = aggregate(PauseScore ~ Codon + Sample, codon_rpm, quantile)
codon_quantiles$PS_min = codon_quantiles$PauseScore[,"0%"]
codon_quantiles$PS_lower = codon_quantiles$PauseScore[,"25%"]
codon_quantiles$PS_median = codon_quantiles$PauseScore[,"50%"]
codon_quantiles$PS_upper = codon_quantiles$PauseScore[,"75%"]
codon_quantiles$PS_max = codon_quantiles$PauseScore[,"100%"]
codon_quantiles = codon_quantiles[,c(1:2,4:ncol(codon_quantiles))]

### SAVE DATA ##################################################################

# Save the table of medians
write.table(
  codon_quantiles,
  "results/2018-04-04/codon_PauseScore_quantiles.with_plasmids.1.tab",
  quote=F, sep="\t", row.names=F, col.names=T
  )
