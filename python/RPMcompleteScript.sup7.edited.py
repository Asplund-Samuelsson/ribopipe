#!/usr/bin/python2.7

from Bio import SeqIO

"""
Supplementary Note 7: Complete RPM-normalized read densities

Author: Annemarie Becker

Modified by: Johannes Asplund-Samuelsson (KTH)

inputFileP:
RPM-normalized read density file for plus strand (Supplementary Note 6)
    col0/NA: reference sequence
    col1/0: position along genome
    col2/1: RPM-normalized read density at that position

inputFileM:
RPM-normalized read density file for minus strand (Supplementary Note 6)
    col0/NA: reference sequence
    col1/0: position along genome
    col2/1: RPM-normalized read density at that position

inputFileG:
Concatenated FASTA file of the reference genome

outputFileP:
complete RPM-normalized read density file for all positions along the genome on plus strand
    col0: reference sequence
    col1: position along genome
    col2: RPM-normalized read density at that position

outputFileM:
complete RPM-normalized read density file for all positions along the genome on minus strand
    col0: reference sequence
    col1: position along genome
    col2: RPM-normalized read density at that position

"""


def RPMcomplete(inputFileP, inputFileM, inputFileG, outputFileP, outputFileM):

    # Initialize RPM dictionaries
    DictP = {}
    DictM = {}

    # Load reference sequence lengths
    ref_lengths = {}

    FastaFile = open(inputFileG, 'rU')

    for rec in SeqIO.parse(FastaFile, 'fasta'):
        # Store length of reference sequence
        ref_lengths[rec.id] = len(rec)
        # Initialize reference sequence entries in RPM dictionaries
        DictP[rec.id] = {}
        DictM[rec.id] = {}

    FastaFile.close()


### PLUS STRAND###

    inFile = open(inputFileP, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        try:
            col0 = str(fields[-3])
        except IndexError:
            # If the reference is not specified, assume fasta has one reference
            col0 = list(DictP.keys())[0]
        col1 = int(fields[-2])
        col2 = float(fields[-1])
        DictP[col0][col1] = col2
        line = inFile.readline()

    outFile = open(outputFileP, 'w')

    for ref in DictP:
        for elem in range(1, ref_lengths[ref] + 1):
            if elem not in DictP[ref]:
                DictP[ref][elem] = 0.0

        ListP = DictP[ref].items()
        ListP.sort()

        for J in ListP:
            outFile.write('\t'.join([ref, str(J[0]), str(J[1])]) + '\n')

    outFile.close()


### MINUS STRAND###

    inFile = open(inputFileM, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split()
        try:
            col0 = str(fields[-3])
        except IndexError:
            col0 = list(DictM.keys())[0]
        col1 = int(fields[-2])
        col2 = float(fields[-1])
        DictM[col0][col1] = col2
        line = inFile.readline()

    outFile = open(outputFileM, 'w')

    for ref in DictM:
        for elem in range(1, ref_lengths[ref] + 1):
            if elem not in DictM[ref]:
                DictM[ref][elem] = 0.0

        ListM = DictM[ref].items()
        ListM.sort()

        for J in ListM:
            outFile.write('\t'.join([ref, str(J[0]), str(J[1])]) + '\n')

    outFile.close()



if __name__=='__main__':
    # Parse commandline arguments
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--inP', help='Input file P.')
    parser.add_argument('--inM', help='Input file M.')
    parser.add_argument('--inG', help='Genome fasta.')
    parser.add_argument('--outP', help='Output file P.')
    parser.add_argument('--outM', help='Output file M.')

    args = parser.parse_args()

    inputFileP = args.inP
    inputFileM = args.inM
    inputFileG = args.inG
    outputFileP = args.outP
    outputFileM = args.outM

    RPMcomplete(inputFileP, inputFileM, inputFileG, outputFileP, outputFileM)
