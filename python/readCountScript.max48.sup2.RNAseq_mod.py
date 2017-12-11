#!/usr/bin/python2.7

"""
Supplementary Note 2: Center Weighting

Authors: Eugene Oh, Annemarie Becker

Modified for RNAseq by: Johannes Asplund-Samuelsson (KTH)

inputFile:
.map file generated by Bowtie default output.

outputFileP:
read density file for plus strand
    col0: position along genome
    col1: read density at that position

outputFileM:
read density file for minus strand
    col0: position along genome
    col1: read density at that position

"""


def rawdata(inputFile, outputFileP, outputFileM):

    pDict = {}
    mDict = {}

    inFile = open(inputFile, 'r')
    line = inFile.readline()
    while line != '':
        fields = line.split("\t")
        col2 = str(fields[1])   #strand; note: if sequencing was performed without barcode reading, the column numbering is changed
        col4 = int(fields[3])   #left-most position
        col5 = str(fields[4])   #read seq
        length = len(col5)      #read length

        if col2 == '+':	#for plus strand
            end5 = col4 + 1 #Bowtie uses zero-based offset, transform to 1-based
            end3 = end5 + length - 1

            for elem in range(end5, end3 + 1):
                if elem in pDict:
                    pDict[elem] += (1.0 / length)
                else:
                    pDict[elem] = (1.0 / length)

        elif col2 == '-': #for minus strand
            end3 = col4 + 1 #for minus strand, Bowtie gives leftmost position (3' end) with zero-based numbering
            end5 = end3 + length - 1

            for elem in range(end3, end5 + 1):
                if elem in mDict:
                    mDict[elem] += (1.0 / length)
                else:
                    mDict[elem] = (1.0 / length)

        line = inFile.readline()

    pList = pDict.items()
    pList.sort()
    outFileP = open(outputFileP, 'w')
    for J in pList:
        outFileP.write(str(J[0]) + '\t' + str(J[1]) + '\n')

    mList = mDict.items()
    mList.sort()
    outFileM = open(outputFileM, 'w')
    for J in mList:
        outFileM.write(str(J[0]) + '\t' + str(J[1]) + '\n')


if __name__=='__main__':
    # Parse commandline arguments
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--infile', help='Input file (bwt).')
    parser.add_argument('--outP', help='Output file P.')
    parser.add_argument('--outM', help='Output file M.')
    parser.add_argument('--min', help='Not used.', type=int) # DUMMY - NOT USED
    parser.add_argument('--max', help='Not used.', type=int) # DUMMY - NOT USED

    args = parser.parse_args()

    inputFile = args.infile
    outputFileP = args.outP
    outputFileM = args.outM

    rawdata(inputFile, outputFileP, outputFileM)
