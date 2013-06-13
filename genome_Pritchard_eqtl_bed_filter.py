#!/usr/bin/env python
import sys
import subprocess

if(len(sys.argv)!= 3):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: Pritchard_eqtl_bed_filter.py [BED In] [BED Out]\n\n')
	sys.exit(1)

filein = open(sys.argv[1],'r')
bedout = open(sys.argv[2],'w')

def name_cleaner( bedline,index ):
	'''Cleans redundancy in bed names.'''
	names = bedline[index].split(';')
	names = set(names)
	return(";".join(names))

def bookend_checker( bedline ):
	'''Fixes lines where mergeBed combines bookended SNPs.'''
	size = int(bedline[2]) - int(bedline[1])
	snps = len(bedline[3].split(';'))
	if size == 1:
		return(0)
	if size > 1 and size == snps:
		return(size)
	if size > 2 and size != snps:
		return(1)
	if size == 2 and snps == 1:
		return('A')

wrongs = 0
wronglines = []
for line in filein:
	liner = line.strip().split('\t')
	eqtl = name_cleaner(liner,3)
	liner[3] = eqtl
	records = bookend_checker(liner)
	if records == 0:
		print >> bedout, "\t".join(liner)
	if records == 1:
		wrongs += 1
		wronglines.append(liner)
		continue
	if records > 1 and records != 'A':
		for time in range(records):
			currline = liner
			currline[1] = str(int(liner[1]) + time)
			currline[2] = str(int(currline[1]) + 1)
			print >> bedout, "\t".join(liner)
	if records == 'A':
		liner[1] = str(int(liner[2]) - 1)
		print >> bedout, "\t".join(liner)

filein.close()
bedout.close()
if wrongs > 0:
	print 'Warning! Will Robinson! There were ' + str(wrongs) + ' error(s)!!'
	if wrongs < 20:
		print wronglines
