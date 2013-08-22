#!/usr/bin/env python
import sys
import subprocess
import gzip

outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/'
combofile = 'EncodeChIP_Uniform_Combined.bed'
counterfile = 'binding_counts.txt'
filekey = open('/mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/Downloads/Uniform_Calls/Target_list.txt','r')
combobed = open(outdir + combofile,'w')
countering = open(outdir + counterfile,'w')
filedict = {}
for line in filekey:
	filedict[line.strip().split()[1]] = line.strip().split()[0]

filekey.close()

for factor in filedict.keys():
	print factor
	counts = 0
	bed = gzip.open('/mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/Downloads/Uniform_Calls/' + filedict[factor],'r')
	for line in bed:
		pos = line.strip().split()[0:3]
		scorepos = [line.strip().split()[6]]
		scorepos.extend(line.strip().split()[8:10])
		print >> combobed, "\t".join(pos) + '\t' + factor + '\t' + "\t".join(scorepos)
		counts += 1
	print counts
	print >> countering, factor + '\t' + str(counts)

combobed.close()
countering.close()

caller = "sort -k1,1 -k2,2n -k3,3n " + outdir + combofile + " > /mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/sorted_" + combofile + "; sort -k1,1 " + outdir + counterfile + " > /mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/sorted_" + counterfile + "; rm " + outdir + combofile + "; rm " + outdir + counterfile
calling = subprocess.Popen(caller,shell=True)
calling.wait()
#remover = "rm " + outdir + combofile
#subprocess.call([remover],shell=True)