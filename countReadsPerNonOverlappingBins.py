#! /usr/bin/python

import sys
import pysam
import matplotlib
import datetime
import numpy as np
import argparse
import re
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--bam',dest='input_bam_file',required=True, help='Input bam file ; required')
parser.add_argument('--binSize',dest='binSize',required=True,help='Size of none overlapping bin in bases; required')
parser.add_argument('--outFile',dest='outputFileName',required=True, help='Name of the output file; required')
parser.add_argument('--RPM',dest='readsPerMillion',type=bool,default=False, help=' Reads Per Million normalization, not default, set it to True if required')
parser.add_argument('--makeplot',dest='makeLinePlot',type=bool,default=False, help='Make line plots for each bin genome wide, not default, set it to True if required')
parser.add_argument('--readSize',dest='manualReadSizeInput',type=bool,default=False, help='Read Size to be input maulally, not default, supply the number if required.\nBy default, the program with take the first 10000 reads length and use the median read length for normalization')

if len(sys.argv)==1:
	parser.print_help(sys.stderr)
	sys.exit(1)

information = parser.parse_args()

binSize = int(information.binSize)
bamFileName = information.input_bam_file
outputFileName = information.outputFileName 
rpmStatus = information.readsPerMillion 
makePlot= information.makeLinePlot
readSizeChecker = information.manualReadSizeInput

def getReadLength(bamfile):
	readSize = []
	flag = 0
	for reads in bamfile:
		if flag >10000:
			break
		else:
			readSize.append(len(reads.seq))
		flag = flag + 1	
	readSize = np.asarray(readSize)
	return np.median(readSize)
	
def getTotalMappedReadsFromBamIndex(bamfileLocation):
	tmp = pysam.idxstats(bamfileLocation)
	mappedReads = []
	for row in tmp[0:-1].split('\n'):
		mappedReads.append(float(row.split('\t')[2]))
	return sum(mappedReads)	

def sorted_nicely( l ):
	""" Sorts the given iterable in the way that is expected.
		https://arcpy.wordpress.com/2012/05/11/sorting-alphanumeric-strings-in-python/
	    Required arguments:
	        l -- The iterable to be sorted."""
		
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key = alphanum_key)

def makeChromosomeIntegerPair(l):
	flag = 0 
	integerChromPair = {}
	for items in l:
		integerChromPair[items] = flag
		flag = flag+1
	return integerChromPair	

def getChromosomeFromBamFile(bamfile):
	chrom = {}
	for lines in bamfile.header['SQ']:
		chrom[lines['SN']] = int(lines['LN'])
	return chrom	


print "Reading bamfile : "+str(bamFileName)
starttime=datetime.datetime.now()

print "starting time : "+str(starttime) 
bamFile=pysam.AlignmentFile(bamFileName,'rb')



chrom = getChromosomeFromBamFile(bamFile)
ChromArray = chrom.keys()

if readSizeChecker == False:
	readLength = getReadLength(bamFile)
else:
	readLength = readSizeChecker

totalMappedReads = getTotalMappedReadsFromBamIndex(bamFileName)
refIntChrom = makeChromosomeIntegerPair(sorted_nicely(ChromArray))


chromRef = []
position = []
#
for read in bamFile:
	chromRef.append(read.reference_id)
	position.append(read.reference_start)	

chromRef = np.asarray(chromRef)
position = np.asarray(position)
print "\n\nReading bamfile complete\nNow Binning the data"
starttime=datetime.datetime.now()
print "starting time : "+str(starttime) +'\n\n'




readCounts = []
chromosomeboundries = [0]	
filename=str(outputFileName)+'.bdg'
bedoutfile=open(filename,'w')
for chromosome in sorted_nicely(ChromArray):
	step=0 # Initial start position
	chromIndexes = np.where(chromRef==refIntChrom[chromosome])
	chromLocation = np.sort(position[chromIndexes], kind='quicksort')
	chromosomeboundries.append((chrom[chromosome]/binSize)+chromosomeboundries[-1])
	step = 0
	flag = 0
	
	while(step <(chrom[chromosome])):
		end=int(step)+int(binSize) #End Coordinate
	
		readCount = np.where(np.logical_and(chromLocation[flag:] >= step,chromLocation[flag:] <= end))
		checkLengthArray = len(readCount[0])
		if end < chrom[chromosome]:
			if rpmStatus == True:
				bedoutfile.write(str(chromosome)+"\t"+str(step)+'\t'+str(end)+'\t'+str((float(len(readCount[0]))*1000000.0)/(totalMappedReads*readLength))+'\n')
				readCounts.append((float(len(readCount[0]))*1000000.0)/(totalMappedReads*readLength))
			else:	
				bedoutfile.write(str(chromosome)+"\t"+str(step)+'\t'+str(end)+'\t'+str(len(readCount[0]))+'\n')
				readCounts.append((float(len(readCount[0]))*1000000.0)/(totalMappedReads*readLength))
		else:
			if rpmStatus == True:
				bedoutfile.write(str(chromosome)+"\t"+str(step)+'\t'+str(chrom[chromosome])+'\t'+str((float(len(readCount[0]))*1000000.0)/(totalMappedReads*readLength))+'\n')
				readCounts.append((float(len(readCount[0]))*1000000.0)/(totalMappedReads*readLength))
			else:
				bedoutfile.write(str(chromosome)+"\t"+str(step)+'\t'+str(chrom[chromosome])+'\t'+str(len(readCount[0]))+'\n')
				readCounts.append((float(len(readCount[0]))*1000000.0)/(totalMappedReads*readLength))
		
	
		if checkLengthArray == 0:
			flag = flag
		else:	
			flag = readCount[0][-1]	
		#print flag	
		step=end
#	print str(chromosome)+" is done"	
bedoutfile.close()  # close file handler
endtime=datetime.datetime.now()

if makePlot == True:
	filename=str(outputFileName)+'.png'
	plt.figure(figsize=(20,15))
	plt.plot(readCounts,color='black')
	maxylim=max(readCounts)+1
	labelLocal = []
	for index in range(len(chromosomeboundries)-1):
		plt.axvline(x=chromosomeboundries[index+1], color='r', linestyle='--')
		labelLocal.append((chromosomeboundries[index]+chromosomeboundries[index+1])/2)	

		#plt.text(((chromosomeboundries[index]+chromosomeboundries[index+1])/2),,sorted_nicely(ChromArray)[index])#,rotation=90)
	plt.xticks(labelLocal,sorted_nicely(ChromArray),rotation=90)	
	plt.axhline(np.median(np.asarray(readCounts)),color='b',linewidth=1)
	plt.suptitle('Whole genome plot\nBin size : '+str(binSize)+'\nGreen line represents the median read distribution', fontsize=16)
	#plt.show()

	plt.savefig(filename,format='png', dpi=1200)
	print "Line plot was saved at : "+str(filename)

print "Process completed at : "+str(endtime)
print "bedgraph file was saved at : "+str(outputFileName
)
