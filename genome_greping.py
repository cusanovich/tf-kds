#!/usr/bin/env python
import os
import sys
import gzip
import subprocess

#if(len(sys.argv)!= 4):
#	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: gff2bed.py [Outdir] [In File] [Out File]\n\n')
#	sys.exit(1)

#motif = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/all_motifs.txt','r')
#motifs = motif.readline().strip().split('|')
genelists = os.listdir('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/')
genelists = sorted(set([x.split('.')[0] for x in genelists]))
factors = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt'


for genes in genelists:
	currtffile = open('../Grepers/'+gene+'.BindingTF.txt','r')
	currdownfile = open('../Grepers/'+gene+'.DownstreamTFs.txt','r')
	currtf = currtffile.readline().strip().split('|')
	currdownstream = currdownstream.readline().strip().split('|')
	currf = open(factors,'r')
	currfactor = []
	currdowns = []
	for line in currf:
		if line.strip().split()[0] in currtf:
			currfactor.append(line.strip().split()[1])
		if line.strip().split()[0] in currdownstream:
			currdowns.append(line.strip().split()[1])

	currf.close()
	currfactor = set(currfactor)
	currdowns = set(currdowns)





gffin = gzip.open(sys.argv[1] + sys.argv[2],'rb')
bedder = sys.argv[1] + sys.argv[3]
bedout = open(bedder,'w')
#if '.bed.gz' not in sys.argv[2]:
#    bedout = gzip.open(sys.argv[2] + '.gz','w')
#else:
#    bedout = gzip.open(sys.argv[2],'w')
chosen = open('../qtlsources.txt','r')

choose = chosen.readlines()
choose = [x.strip() for x in choose]
chosen.close()

for line in gffin:
    liner = line.strip().split()
    if liner[1] not in choose:
        continue
    chrom = liner[0]
    start = str(int(liner[3]) - 1)
    end = liner[4]
    name = liner[2]
    score = liner[5]
    print >> bedout, chrom + '\t' + start + '\t' + end + '\t' + name + '\t' + score

gffin.close()
bedout.close()

outdir = sys.argv[1]
sorter = 'sort -k1,1 -k2,2n -k3,3n -k4,4 ' + bedder + ' > ' + outdir + 'sorted_' + sys.argv[3] + '; gzip ' + outdir + 'sorted_' + sys.argv[3]
subprocess.Popen(sorter,shell=True)
