#!/usr/bin/python2.7

"""
Supplementary Note 6: RPM-normalized read densities

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

inputNumber:
total read number as float (Supplementary Note 3)

outputFileP:
RPM-normalized read density file for plus strand
    col0/NA: reference sequence
    col1/0: position along genome
    col2/1: RPM-normalized read density at that position

outputFileM:
RPM-normalized read density file for minus strand
    col0/NA: reference sequence
    col1/0: position along genome
    col2/1: RPM-normalized read density at that position

"""


def norm(inputFileP, inputFileM, inputNumber, outputFileP, outputFileM):

### PLUS STRAND ###

    inFile = open(inputFileP, 'r')
    inNumber = open(inputNumber, 'r')
    outFile = open(outputFileP, 'w')

    line = inFile.readline()
    number = inNumber.readline()
    totalReads = int(float(number))

    i = 0
    while line != '':
        fields = line.split()
        rc = float(fields[-1])
        RPM = rc / totalReads * 1000000
        outFile.write('\t'.join(fields[0:-1]) + '\t' + str(RPM) + '\n')
        line = inFile.readline()


### MINUS STRAND ###

    inFile = open(inputFileM, 'r')
    inNumber = open(inputNumber, 'r')
    outFile = open(outputFileM, 'w')

    line = inFile.readline()
    number = inNumber.readline()
    totalReads = int(float(number))

    i = 0
    while line != '':
        fields = line.split()
        rc = float(fields[-1])
        RPM = rc / totalReads * 1000000
        outFile.write('\t'.join(fields[0:-1]) + '\t' + str(RPM) + '\n')
        line = inFile.readline()



if __name__=='__main__':
    # Parse commandline arguments
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--inP', help='Input file P.')
    parser.add_argument('--inM', help='Input file M.')
    parser.add_argument('--number', help='Input number.')
    parser.add_argument('--outP', help='Output file P.')
    parser.add_argument('--outM', help='Output file M.')

    args = parser.parse_args()

    inputFileP = args.inP
    inputFileM = args.inM
    inputNumber = args.number
    outputFileP = args.outP
    outputFileM = args.outM

    norm(inputFileP, inputFileM, inputNumber, outputFileP, outputFileM)
