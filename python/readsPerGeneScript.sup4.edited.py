#!/usr/bin/python2.7

from Bio import SeqIO

"""
Supplementary Note 4: Read density per gene

Authors: Eugene Oh

Modified by: Johannes Asplund-Samuelsson (KTH)

inputFileP:
read density file for plus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

inputFileM:
read density file for minus strand (Supplementary Note 2)
    col0: position along genome
    col1: read density at that position

inputListP:
E. coli MC4100 gene list of the plus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene

inputListM
E. coli MC4100 gene list of the minus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene

outputFileP:
read densities per gene on plus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: sum of read densities

outputFileM:
read densities per gene on minus strand
    col0: gene name
    col1: start coordinate of gene
    col2: stop coordinate of gene
    col3: sum of read densities

genomeLength:
length of genome sequence for negative gene position handling
    int

"""


def expression(inputFileP, inputFileM, inputListP, inputListM, outputFileP, \
outputFileM, inputFileG):

    # Initialize readcount dictionaries
    DictP = {}
    DictM = {}

    # Load reference sequence lengths
    ref_lengths = {}

    FastaFile = open(inputFileG, 'rU')

    for rec in SeqIO.parse(FastaFile, 'fasta'):
        # Store length of reference sequence
        ref_lengths[rec.id] = len(rec)
        # Initialize reference sequence entries in readcount dictionaries
        DictP[rec.id] = {}
        DictM[rec.id] = {}

    FastaFile.close()

### PLUS STRAND ###

  # Upload read density file from plus strand as a dictionary
    def load_input(infile, Dict):
        inFile = open(infile, 'r')
        line = inFile.readline()
        while line != '':
            fields = line.split()
            try:
                col0 = str(fields[-3])
            except IndexError:
                # If no reference, assume fasta has one reference
                col0 = list(Dict.keys())[0]
            col1 = int(fields[-2])
            col2 = float(fields[-1])
            Dict[col0][col1] = col2
            line = inFile.readline()

    load_input(inputFileP, DictP)

  # Upload plus strand gene list as a dictionary and list

    def assign_gene_reads(inputList, Dict, outputFile):

        geneDict = {} # create dictionary with col0=gene name; col1=read number
        geneList = [] # create list that looks like input gene list

        inFile = open(inputList, 'r')
        line = inFile.readline()
        while line != '':
            fields = line.split()
            geneList.append(fields)        # add an item to the end of the list
            ref = str(fields[0])            # reference sequence ID
            gene = str(fields[1])           # gene name
            start = int(fields[2])          # start
            stop = int(fields[3])           # stop

      # Sum up and write read densities per protein coding region in dictionary

            for Z in range(start, stop + 1):
                if Z < 0:
                    # Handle gene positions that are negative or zero
                    Z = ref_lengths[ref] + Z + 1
                if Z in Dict[ref] and gene in geneDict:
                    geneDict[gene] += Dict[ref][Z]
                elif Z in Dict[ref]:
                    geneDict[gene] = Dict[ref][Z]
            line = inFile.readline()

      # Assign gene expression levels to all genes

        tupledlist = geneDict.items()
        for J in geneList:
            match = 0
            for K in tupledlist:
                if J[1] == K[0]:
                    match = 1
                    J.append(K[1])
            if match == 0:		#list genes that don't have any reads
                J.append(0)

      # Output file for plus strand

        outFile = open(outputFile, 'w')
        for J in geneList:
            output = '\t'.join([str(x) for x in J]) + '\n'
            outFile.write(output)

    assign_gene_reads(inputListP, DictP, outputFileP)

### MINUS STRAND ###

  # Upload read density file from minus strand as a dictionary

    load_input(inputFileM, DictM)

  # Upload minus strand gene list as a dictionary and list

    assign_gene_reads(inputListM, DictM, outputFileM)


if __name__ == '__main__':
    # Parse commandline arguments
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--inP', help='Input file P.')
    parser.add_argument('--inM', help='Input file M.')
    parser.add_argument('--listP', help='Input list P.')
    parser.add_argument('--listM', help='Input list M.')
    parser.add_argument('--outP', help='Output file P.')
    parser.add_argument('--outM', help='Output file M.')
    parser.add_argument('--inG', help='Genome fasta.')

    args = parser.parse_args()

    inputFileP = args.inP
    inputFileM = args.inM
    inputListP = args.listP
    inputListM = args.listM
    outputFileP = args.outP
    outputFileM = args.outM
    inputFileG = args.inG

    expression(inputFileP, inputFileM, inputListP, inputListM, outputFileP, outputFileM, inputFileG)
