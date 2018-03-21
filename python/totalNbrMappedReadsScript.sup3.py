#!/usr/bin/python2.7

"""
Supplementary Note 3: Total reads

Author: Annemarie Becker

Modified by: Johannes Asplund-Samuelsson (KTH)

inputFileP:
read density file for plus strand (Supplementary Note 2)
    col0/NA: reference sequence
    col1/0: position along genome
    col2/1: read density at that position

inputFileM:
read density file for minus strand (Supplementary Note 2)
    col0/NA: reference sequence
    col1/0: position along genome
    col2/1: read density at that position

outputFile:
total read number as float

"""


def countReads(inputFileP,inputFileM, outputFile):

    inFileP = open(inputFileP, 'r')
    inFileM = open(inputFileM, 'r')
    outFile = open(outputFile, 'w')

    line = inFileP.readline()
    i = 0
    while line != '':
        fields = line.split()
        rc = float(fields[-1])
        i = i + rc # count reads on plus strand
        line = inFileP.readline()
    totReadsP = i

    line = inFileM.readline()
    j = 0
    while line != '':
        fields = line.split()
        rc = float(fields[-1])
        j = j + rc # count reads on minus strand
        line = inFileM.readline()
    totReadsM = j

    totalReads = i + j
    outFile.write(str(totalReads))


if __name__ == '__main__':
    # Parse commandline arguments
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--inP', help='Input file P.')
    parser.add_argument('--inM', help='Input file M.')
    parser.add_argument('--out', help='Output file.')

    args = parser.parse_args()

    inputFileP = args.inP
    inputFileM = args.inM
    outputFile = args.out

    countReads(inputFileP,inputFileM, outputFile)
