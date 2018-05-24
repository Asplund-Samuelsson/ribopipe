#!/usr/bin/env python3

# Import libraries
import re

# Define infiles
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--inGen', help='Genome fasta.')
parser.add_argument('--inPos', help='Positions file.')

args = parser.parse_args()

genome_fasta = args.inGen
pos_file = args.inPos

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

# Load the positions for each 5' UTR
all_utr_positions = {}

for line in open(pos_file).readlines():
    line = line.strip().split("\t")
    if line[0] == "Name":
        continue
    utr_id = "\t".join(line[0:3])
    position = int(line[3]) - 1 # Python starts from 0, genome starts from 1
    try:
        all_utr_positions[utr_id].append(position)
    except KeyError:
        all_utr_positions[utr_id] = [position]

# Go through each UTR and extract its sequence
for utr_id in all_utr_positions:
    # Get ORF name, sequence and strand
    orf, seq, strand = utr_id.split("\t")
    # Get nucleotide positions in order
    positions = sorted(all_utr_positions[utr_id])
    # Get the nucleotide sequence
    UTR = "".join([genome_sequence[seq][x] for x in positions])
    # Turn sequence into reverse complement if it is on the minus strand
    if strand == "-":
        UTR = revcomp(UTR)
    # Print in FASTA format
    print(">" + orf)
    print(UTR)
