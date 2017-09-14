#!/usr/bin/python3.5

#define inputFile from commandline argument
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', help='inputFile.')

args = parser.parse_args()

inputFile = args.i
dotsplit = inputFile.split('.')
outputFile = dotsplit[0] + "." + dotsplit[1] + ".lengthDist"

inFile = open(inputFile, 'r') #create a handle (link called "outFile") to a readable ('r') file
Dict = {}
line = inFile.readline()

while line != '':			#read row in text file
	fields = line.split()	#create and array with fields in the row as elements
	length = len(fields[5]) #read length
	if length in Dict:		#if item "length" exist in Dict
		Dict[length] += 1	#add 1
	else:					#otherwise
		Dict[length] = 1	#add new item with value 1
	line = inFile.readline()

List = Dict.items()
outFile = open(outputFile, 'w')	#write one row for each item in List: item[0] = first tuple, (length), item[1] = second tuple (count), '\n' = new row
outFile.write("read length (bp)" + '\t' + "count" + '\n')	#write header, '\t' = tab

for item in List:
	outFile.write(str(item[0]) + '\t' + str(item[1]) + '\n')	#write one row for each item in List: item[0] = first tuple, (length), item[1] = second tuple (count), '\n' = new row

outFile.close()


#Plot data and save as png
x_list = []
y_list = []
for item in List:
	x_list.append(item[0])	#fill x-list with lengths
	y_list.append(item[1])	#fill y-list with counts

import matplotlib.pyplot as plt
plt.plot(x_list, y_list, 'r-')
plt.grid(True)
plt.savefig(outputFile + ".png")

