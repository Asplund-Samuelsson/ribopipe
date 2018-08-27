#!/usr/bin/env Rscript

### FILENAMES, SAMPLE IDS AND STRAND IDS #######################################

# Define input directory
indir="."
args = commandArgs(trailingOnly=T)
genelist_file = args[1]
shift = args[2]
# indir = "/hdd/common/proj/RibosomeProfiling/results/2018-04-16/CSD2_seqmagick"
# genelist_file = "tools/ribopipe/genelists/syn_PCC6803/NC_000911.1_chr_7plasmids.genelist_full.tab"
# shift = -12

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
gen_tc_data = unlist(lapply(gen_tc_files, scan))

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

# Load gene types
gene_types = read.table(
  genelist_file, stringsAsFactors=F, sep="\t", header=T, quote=""
)[,c("Type","Old_locus_tag", "Sequence")]
colnames(gene_types)[grep("Old", colnames(gene_types))] = "Name"

### SHIFT RPM VALUES ###########################################################

# SIGNAL SHIFT
# shift = -12

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

### RECALCULATE RPM PER GENE ###################################################

# Remove current read count
gen = gen[,grep("(Sample)|(Reads)", colnames(gen), invert=T, perl=T)]
gen = unique(gen)

# Subset to coding sequences only
gen = subset(gen, Name %in% subset(gene_types, Type == "CDS")$Name)

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

# Calculate gene length
gen$Length = gen$End - gen$Start + 1

# Special treatment for genes with negative start
gen$Length = ifelse(gen$Start > 0, gen$Length, gen$End - gen$Start)

# Calculate RPKM
gen$RPKM = 1000 * gen$RPM / gen$Length

# Reconstruct per-position read count
rpm = inner_join(rpm, tc_gen)
rpm$Reads = rpm$RPM * rpm$TotalReads / 1000000

# Calculate RPK
gen$RPK_gene = 1000 * gen$Reads / gen$Length
rpm$RPK_nucl = 1000 * rpm$Reads

### ANALYZE ####################################################################

# Expand each gene to all positions
gen_allpos = as.data.frame(rbindlist(lapply(
  unique(gen$Name),
  function (x) {
    gen_s = subset(gen, Name == x)
    Position = (unique(gen_s$Start) - 50):(unique(gen_s$End) + 50)
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
gen_allpos = inner_join(
  gen_allpos, rpm[,c("Sequence","strand","Position","Sample","RPK_nucl")]
  )

# Calculate PauseScore
gen_allpos$PauseScore = gen_allpos$RPK_nucl / gen_allpos$RPK_gene

# Calculate relative position
gen_allpos$RelPos = ifelse(
  gen_allpos$strand == "+",
  gen_allpos$Position - gen_allpos$Start,
  gen_allpos$End - gen_allpos$Position
  )

# Calculate relative position for genes crossing position 1
gen_allpos$RelPos = ifelse(
  gen_allpos$Start < 1,
  # Treat only genes crossing position 1
  ifelse(
      gen_allpos$Position > gen_allpos$End + 50,
      # Return to negative positions if before 1
      ifelse(
        gen_allpos$strand == "+",
        (gen_allpos$Position - gen_allpos$Size) - gen_allpos$Start - 1,
        gen_allpos$End - (gen_allpos$Position - gen_allpos$Size)
        ),
      ifelse(
        gen_allpos$strand == "+",
        # Plus-strand genes should end on rel. pos. length - 1 (starts at 0)
        gen_allpos$Position - gen_allpos$Start - 1,
        gen_allpos$End - gen_allpos$Position
        )
    ),
  # If the gene does not cross position 1, keep relative position as it is
  gen_allpos$RelPos
  )

# Rename to gene_profiles
gene_profiles = gen_allpos

# Calculate distance from end of gene
gene_profiles$RelPosFromEnd = gene_profiles$RelPos - gene_profiles$Length + 1

# Find positions with multiple ORFs
orf_positions = unique(
  gene_profiles[,c("Sequence", "strand", "Name", "Position")]
  )

multi_orf_positions = orf_positions[
  duplicated(orf_positions[,c("Sequence", "strand", "Position")]),
  c("Sequence", "strand", "Position")
]

multi_orf_positions_with_RelPos = inner_join(
  multi_orf_positions,
  unique(gene_profiles[,
    c("Sequence", "strand", "Name", "Position", "RelPos", "RelPosFromEnd")
    ])
  )

# Calculate whether the position is within the ORF or not
multi_orf_positions_with_RelPos$WithinORF = multi_orf_positions_with_RelPos$RelPos >= 0 & multi_orf_positions_with_RelPos$RelPosFromEnd <= 0

# Determine which positions should be banned for each ORF
positions_within_orf = filter(
  multi_orf_positions_with_RelPos, WithinORF
)[,c("Sequence", "strand", "Name", "Position")]

colnames(positions_within_orf)[3] = "OtherName"

multi_orf_positions_with_RelPos = inner_join(
  multi_orf_positions_with_RelPos, positions_within_orf
)

orf_positions_that_overlap_with_other_orfs = filter(
  multi_orf_positions_with_RelPos, Name != OtherName
)

bad_positions = orf_positions_that_overlap_with_other_orfs[,
  c("Sequence", "strand", "Name", "Position")
]

# Subset gene profiles to positions that are not in multiple ORFs
gene_profiles = anti_join(gene_profiles, bad_positions)

# Save results as table
write.table(
  gene_profiles,
  "analysis/gene_PauseScore_profiles.all_genes.tab",
  quote=F, sep="\t", row.names=F, col.names=T
)

# Remove gene that skews 5' UTR
gene_profiles = filter(gene_profiles, Name != "slr1474")

# Remove genes with 128 reads or less
genes_below_threshold = as.character(unique(filter(gene_profiles, Reads < 128)$Name))
gene_profiles = subset(gene_profiles, !(Name %in% genes_below_threshold))

# Subset for 5 and 3 prime ends
gene_prof_5prime = filter(gene_profiles, RelPos %in% -50:150)
gene_prof_3prime = filter(gene_profiles, RelPosFromEnd %in% -150:50)

# Calculate means
mean_profile_5prime = aggregate(PauseScore ~ Sample + RelPos, gene_prof_5prime, mean)
mean_profile_3prime = aggregate(PauseScore ~ Sample + RelPosFromEnd, gene_prof_3prime, mean)

### PLOT #######################################################################

library(ggplot2)

# Combine into one dataframe
mean_profile_3prime$Part = "3_prime"
mean_profile_5prime$Part = "5_prime"

colnames(mean_profile_3prime)[2] = "RelPos"

average_profiles = rbind(mean_profile_5prime, mean_profile_3prime)

average_profiles$Part = factor(average_profiles$Part, levels = c("5_prime","3_prime"))

# Mean lines of all genes, by part of ORF, faceted by sample
gp = ggplot(average_profiles, aes(x=RelPos, y=PauseScore))
gp = gp + geom_line(size=0.3)
gp = gp + theme_bw()
gp = gp + facet_grid(Sample ~ Part, scales="free_x")
gp = gp + scale_y_log10()
gp = gp + scale_x_continuous(minor_breaks=seq(-174, 168, 3))
gp = gp + theme(strip.background = element_blank())

n_samples = length(unique(average_profiles$Sample))

ggsave("analysis/average_gene_profile.per_sample.pdf", gp, height=6.875*n_samples*0.25, width=6.875)
# ggsave("art/2018-08-27/average_gene_profile.per_sample.test_1.pdf", gp, height=6.875*n_samples*0.25, width=6.875)
