#!/usr/bin/env python
import os
import sys
import random
import numpy
from kdfunc import genereader, pathmaker

windowsize = '10kb'
grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/Union/' + windowsize + '/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/Union/' + windowsize + '/'
permdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Perms/Union/' + windowsize + '/'
genelists = os.listdir(grepdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))
#genelist = ['ARNTL2_old']
factors = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt'
masterfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowsize + '_combined_midpoint_new_sorted.bed'
centifile = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'results_phastcons_sorted_midpoint.bed'

for gene in genelist:
#for gene in ['YY1_second']:
	print gene
	sys.stdout.flush()
	currdownstream = genereader(grepdir+gene+'.DownstreamTFs.txt')
	currboundplus = genereader(grepdir+gene+'.BoundPlus.txt')
	currdeplus = genereader(grepdir+gene+'.DEandBoundPlus.txt')
	currf = open(factors,'r')
	currdowns = []
	for line in currf:
		if line.strip().split()[0] in currdownstream:
			currdowns.append(line.strip().split()[0])
	currf.close()
	currdowns = set(currdowns)
	print 'Splitting Master Records...'
	sys.stdout.flush()
	master = open(masterfile,'r')
	if len(currdowns) == 0:
		continue
	tenperm = open(outdir + 'GrepDistance/' + gene + '.down.permadistance','w')
	for line in master:
		strand = 1
		liner = line.strip().split()
		if len(liner) < 7:
			continue
		if liner[5] == '-':
			strand = -1
		if liner[9] in currdowns:
			if liner[3] in currboundplus or liner[3] in currdeplus:
				print >> tenperm, liner[3] + "\t" + str(strand*(int(liner[8]) - int(liner[2])))
	master.close()
	tenperm.close()