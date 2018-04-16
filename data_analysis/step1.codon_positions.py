#!/usr/bin/env python3

# Import libraries
import re

# Define infiles
genome_fasta = "/ssd/common/tools/ribopipe/reference_sequences/syn_PCC6803/NC_000911.1_chr_7plasmids.fasta"
CDS_file = "data/2018-03-28/Syn_PCC6803.1_chr_7plasmids.no_pseudo.CDS.tab"

# Load genome sequence
genome_sequence = {}

for line in open(genome_fasta).readlines():
    if line[0] == ">":
        line = line.split()
        sequence = re.sub(">", "", line[0])
        genome_sequence[sequence] = []
        continue
    genome_sequence[sequence].append(line.strip())

for sequence in genome_sequence.keys():
    genome_sequence[sequence] = "".join(genome_sequence[sequence])

# Define reverse complement function
def revcomp(sequence):
    rc_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
    reverse_complement = []
    for nucleotide in sequence[::-1]:
        reverse_complement.append(rc_dict[nucleotide])
    return "".join(reverse_complement)

# Go through each coding sequence, extract codons and positions
for line in open(CDS_file).readlines():
    # Split line and extract components
    line = line.strip().split("\t")
    seq = line[0]
    start = int(line[1])
    end = int(line[2])
    gene = line[3]
    strand = line[4]
    # Get the nucleotide sequence
    if start < 1:
        # If the start is lower than 1, it begins at the end of the genome
        CDS = genome_sequence[seq][start:] + genome_sequence[seq][:end]
    else:
        # Python counts from 0, while first position in genome is designated 1
        CDS = genome_sequence[seq][start-1:end]
    # Set up list of genomic positions for all nucleotides in the sequence
    positions = list(filter(lambda x: x != 0, range(start,end+1)))
    # Modify list
    for i in range(len(positions)):
        if positions[i] < 0:
            positions[i] = len(genome_sequence[seq]) + positions[i] + 1
        else:
            continue
    # Turn sequence into reverse complement if it is on the minus strand
    if strand == "-":
        CDS = revcomp(CDS)
        # ...and reverse the list of positions
        positions.reverse()
    # Go through each codon and print its start and end position
    n = 0
    while n < len(CDS):
        codon = CDS[n:n+3]
        # Print positions of each nucleotide in the codon
        output = "\t".join([
            seq, gene, codon, str(n+1), strand,
            str(positions[n]), str(positions[n + 1]), str(positions[n + 2])
            ])
        print(output)
        n += 3
