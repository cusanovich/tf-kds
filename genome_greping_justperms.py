#!/usr/bin/env python
import os
import sys
import random
import numpy

indir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/'
countfile = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PermCounts.txt','r')
genecounts = {}
for line in countfile:
	liner = line.strip().split()
	genecounts[liner[0]] = [int(x) for x in liner[1:]]

countfile.close()
print 'Making random estimates...'
sys.stdout.flush()
tests = os.listdir(indir)
x=1
for test in sorted(tests):
#for test in ['GrepChipvCenti']:
#	if 'Distance' not in test:
#		continue
	if 'Grep' not in test:
		continue
	masterrandostf = {}
	masterrandosdown = {}
	for i in range(1000):
		masterrandostf[i] = []
		masterrandosdown[i] = []
	permresults = open(outdir + test + '.master.perms','w')
	genes = os.listdir(indir + test)
	for gene in sorted(genes):
#	for gene in ['YY1_second.down.perm.cage']:
#		if 'YY1' not in gene:
#			continue
		if '.perm.' not in gene:
			continue
		if '.tf.' in gene:
			currflen = genecounts[gene.split('.')[0]][0]
			currclen = genecounts[gene.split('.')[0]][2]
			tfscope = 'TF'
		if '.down.' in gene:
			currflen = genecounts[gene.split('.')[0]][1]
			currclen = genecounts[gene.split('.')[0]][3]
			tfscope = 'Downstream'
		justgene = gene.split('.')[0]
		print gene + ' randoms (' + str(currflen) + '/' + str(currclen) + ' genes)...'
		sys.stdout.flush()
		currdict = {}
		currperm = open(indir + test + '/' + gene,'r')
		for line in currperm:
			if line == '':
				print 'pass'
				continue
			liner = line.strip().split()
			if liner[0] not in currdict.keys():
				currdict[liner[0]] = [liner[1]]
			else:
				currdict[liner[0]].append(liner[1])
		currperm.close()
		print 'Current dictionary has ' + str(len(currdict)) + ' genes...'
		sys.stdout.flush()
		#print currdict[currdict.keys()[0]]
		generand = []
		for g in range(1000):
			#print g
			random.seed(x)
		#	if 'IRF4' in gene:
		#		randogenes = sorted(currdict.keys()).extend(random.sample(sorted(currdict.keys()),50))
			if 'perm.p' not in gene:
				randogenes = random.sample(sorted(currdict.keys()),currflen)
			if 'perm.p' in gene:
				randogenes = random.sample(sorted(currdict.keys()),currclen)
			randos = []
			for rando in randogenes:
				randos.extend(currdict[rando])
			if 'GrepChip' in test:
				if '.tf.' in gene:
					masterrandostf[g].extend(randos)
				if '.down.' in gene:
					masterrandosdown[g].extend(randos)
				randy = randos.count('centi')/float(len(randos))
				generand.append(str(randy))
				x += 1
				continue
			if 'GrepDistance' in test:
				randy = [abs(int(z)) for z in randos]
				if '.tf.' in gene:
					masterrandostf[g].extend(randy)
				if '.down.' in gene:
					masterrandosdown[g].extend(randy)
				generand.append(str(numpy.median(randy)))
				x += 1
				continue
			randy = [float(z) for z in randos]
			if '.tf.' in gene:
				masterrandostf[g].extend(randy)
			if '.down.' in gene:
				masterrandosdown[g].extend(randy)
			generand.append(str(numpy.median(randy)))
			x += 1
		print >> permresults, justgene + "\t" + tfscope + "\t" + "\t".join(generand)
	permresults.close()
	currpermresultstf = open(outdir + test + '.alt.tf.perms','w')
	currpermresultsdown = open(outdir + test + '.alt.down.perms','w')
	for rep in range(1000):
		if 'GrepChip' in test:
			print >> currpermresultstf, str(masterrandostf[rep].count('centi')/float(len(masterrandostf[rep])))
			print >> currpermresultsdown, str(masterrandosdown[rep].count('centi')/float(len(masterrandosdown[rep])))
		else:
			print >> currpermresultstf, str(numpy.median(masterrandostf[rep]))
			print >> currpermresultsdown, str(numpy.median(masterrandosdown[rep]))
	currpermresultstf.close()
	currpermresultsdown.close()