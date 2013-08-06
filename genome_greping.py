#!/usr/bin/env python
import os
import sys
import random
import numpy

windowsize = '1kb'
grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/' + windowsize + '/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/' + windowsize + '/'
permdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Perms/' + windowsize + '/'
genelists = os.listdir(grepdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))
#genelist = ['ARNTL2_old']
factors = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt'
masterfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowsize + '_combined_midpoint_sorted.bed'
centifile = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'results_phastcons_sorted_midpoint.bed'

def genereader(thisfile):
	currfile = open(thisfile,'r')
	return(currfile.readline().strip().split('|'))

def pathmaker(path):
	try:
		os.makedirs(path)
	except OSError:
		print 'Warning: ' + path + ' already exists!'

genecounts = {}
pathmaker(outdir + 'GrepDistance/')
pathmaker(outdir + 'GrepCage/')
pathmaker(outdir + 'GrepChipvCenti/')
pathmaker(outdir + 'GrepPhastCons/')
pathmaker(outdir + 'GrepPWM/')
pathmaker(outdir + 'GrepPosterior/')
for gene in genelist:
#for gene in ['YY1_second']:
	print gene
	sys.stdout.flush()
	currtf = genereader(grepdir+gene+'.BindingTF.txt')
	currdownstream = genereader(grepdir+gene+'.DownstreamTFs.txt')
	currbound = genereader(grepdir+gene+'.Bound.txt')
	currboundplus = genereader(grepdir+gene+'.BoundPlus.txt')
	currde = genereader(grepdir+gene+'.DEandBound.txt')
	currdeplus = genereader(grepdir+gene+'.DEandBoundPlus.txt')
	first = 0
	second = 0
	if currde != ['']:
		first = len(currde)
	if currdeplus != ['']:
		second = len(currdeplus)
	genecounts[gene] = [first,second,0,0]
	genecounts[gene] = [len(currde),len(currdeplus),0,0]
	currf = open(factors,'r')
	currfactor = []
	currdowns = []
	currcentifactor = []
	currcentidowns = []
	for line in currf:
		if line.strip().split()[0] in currtf:
			currfactor.append(line.strip().split()[1])
			if len(line.strip().split()[1].split('_')) == 1:
				currcentifactor.append(line.strip().split()[1])
		if line.strip().split()[0] in currdownstream:
			currdowns.append(line.strip().split()[1])
			if len(line.strip().split()[1].split('_')) == 1:
				currcentidowns.append(line.strip().split()[1])
	currf.close()
	currfactor = set(currfactor)
	currdowns = set(currdowns)
	currcentifactor = set(currcentifactor)
	currcentidowns = set(currcentidowns)
	print 'Splitting Master Records...'
	sys.stdout.flush()
	master = open(masterfile,'r')
	tfcagegenes = []
	downcagegenes = []
	if len(currfactor) > 0:
		one = open(outdir + 'GrepDistance/' + gene + '.tf.bound.distance','w')
		oneperm = open(outdir + 'GrepDistance/' + gene + '.tf.perm.distance','w')
		two = open(outdir + 'GrepCage/' + gene + '.tf.bound.cage','w')
		twoperm = open(outdir + 'GrepCage/' + gene + '.tf.perm.cage','w')
		three = open(outdir + 'GrepChipvCenti/' + gene + '.tf.bound.chipvcenti','w')
		threeperm = open(outdir + 'GrepChipvCenti/' + gene + '.tf.perm.chipvcenti','w')
		oneone = open(outdir + 'GrepDistance/' + gene + '.tf.debound.distance','w')
		twotwo = open(outdir + 'GrepCage/' + gene + '.tf.debound.cage','w')
		threethree = open(outdir + 'GrepChipvCenti/' + gene + '.tf.debound.chipvcenti','w')
	if len(currdowns) > 0:
		ten = open(outdir + 'GrepDistance/' + gene + '.down.boundplus.distance','w')
		tenperm = open(outdir + 'GrepDistance/' + gene + '.down.perm.distance','w')
		eleven = open(outdir + 'GrepCage/' + gene + '.down.boundplus.cage','w')
		elevenperm = open(outdir + 'GrepCage/' + gene + '.down.perm.cage','w')
		twelve = open(outdir + 'GrepChipvCenti/' + gene + '.down.boundplus.chipvcenti','w')
		twelveperm = open(outdir + 'GrepChipvCenti/' + gene + '.down.perm.chipvcenti','w')
		tenten = open(outdir + 'GrepDistance/' + gene + '.down.deboundplus.distance','w')
		eleveneleven = open(outdir + 'GrepCage/' + gene + '.down.deboundplus.cage','w')
		twelvetwelve = open(outdir + 'GrepChipvCenti/' + gene + '.down.deboundplus.chipvcenti','w')
	for line in master:
		strand = 1
		liner = line.strip().split()
		if len(liner) < 7:
			continue
		if liner[5] == '-':
			strand = -1
		if liner[9] in currfactor:
			if liner[3] in currbound:
				print >> one, strand*(int(liner[8]) - int(liner[2]))
				print >> oneperm, liner[3] + "\t" + str(strand*(int(liner[8]) - int(liner[2])))
				if liner[3] not in tfcagegenes:
					print >> two, liner[4].split('_')[0]
					print >> twoperm, liner[3] + "\t" + liner[4].split('_')[0]
					tfcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> three, chiper
				print >> threeperm, liner[3] + "\t" + chiper
			if liner[3] in currde:
				print >> oneone, strand*(int(liner[8]) - int(liner[2]))
				if liner[3] not in tfcagegenes:
					print >> twotwo, liner[4].split('_')[0]
					tfcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> threethree, chiper
		if liner[9] in currdowns:
			if liner[3] in currboundplus:
				print >> ten, strand*(int(liner[8]) - int(liner[2]))
				print >> tenperm, liner[3] + "\t" + str(strand*(int(liner[8]) - int(liner[2])))
				if liner[3] not in downcagegenes:
					print >> eleven, liner[4].split('_')[0]
					print >> elevenperm, liner[3] + "\t" + liner[4].split('_')[0]
					downcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> twelve, chiper
				print >> twelveperm, liner[3] + "\t" + chiper
			if liner[3] in currdeplus:
				print >> tenten, strand*(int(liner[8]) - int(liner[2]))
				if liner[3] not in downcagegenes:
					print >> eleveneleven, liner[4].split('_')[0]
					downcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> twelvetwelve, chiper
	master.close()
	if len(currfactor) > 0:
		one.close(),oneperm.close(),two.close(),twoperm.close(),three.close(),threeperm.close(),oneone.close(),twotwo.close(),threethree.close()
	if len(currdowns) > 0:
		ten.close(),tenperm.close(),eleven.close(),elevenperm.close(),twelve.close(),twelveperm.close(),tenten.close(),eleveneleven.close(),twelvetwelve.close()
	print 'Splitting Centipede Records...'
	sys.stdout.flush()
	master = open(centifile,'r')
	if len(currcentifactor) > 0:
		one = open(outdir + 'GrepPhastCons/' + gene + '.tf.bound.phast','w')
		oneperm = open(outdir + 'GrepPhastCons/' + gene + '.tf.perm.phast','w')
		two = open(outdir + 'GrepPWM/' + gene + '.tf.bound.pwm','w')
		twoperm = open(outdir + 'GrepPWM/' + gene + '.tf.perm.pwm','w')
		three = open(outdir + 'GrepPosterior/' + gene + '.tf.bound.post','w')
		threeperm = open(outdir + 'GrepPosterior/' + gene + '.tf.perm.post','w')
		oneone = open(outdir + 'GrepPhastCons/' + gene + '.tf.debound.phast','w')
		twotwo = open(outdir + 'GrepPWM/' + gene + '.tf.debound.pwm','w')
		threethree = open(outdir + 'GrepPosterior/' + gene + '.tf.debound.post','w')
	if len(currcentidowns) > 0:
		ten = open(outdir + 'GrepPhastCons/' + gene + '.down.boundplus.phast','w')
		tenperm = open(outdir + 'GrepPhastCons/' + gene + '.down.perm.phast','w')
		eleven = open(outdir + 'GrepPWM/' + gene + '.down.boundplus.pwm','w')
		elevenperm = open(outdir + 'GrepPWM/' + gene + '.down.perm.pwm','w')
		twelve = open(outdir + 'GrepPosterior/' + gene + '.down.boundplus.post','w')
		twelveperm = open(outdir + 'GrepPosterior/' + gene + '.down.perm.post','w')
		tenten = open(outdir + 'GrepPhastCons/' + gene + '.down.deboundplus.phast','w')
		eleveneleven = open(outdir + 'GrepPWM/' + gene + '.down.deboundplus.pwm','w')
		twelvetwelve = open(outdir + 'GrepPosterior/' + gene + '.down.deboundplus.post','w')
	currcentigenes = []
	currcentidowngenes = []
	for line in master:
		liner = line.strip().split()
		if len(liner) < 7:
			continue
		if liner[9] in currcentifactor:
			if liner[3] in currbound:
				print >> one, liner[13]
				print >> oneperm, liner[3] + "\t" + liner[13]
				print >> two, liner[12]
				print >> twoperm, liner[3] + "\t" + liner[12]
				print >> three, liner[10]
				print >> threeperm, liner[3] + "\t" + liner[10]
			if liner[3] in currde:
				currcentigenes.append(liner[3])
				print >> oneone, liner[13]
				print >> twotwo, liner[12]
				print >> threethree, liner[10]
		if liner[9] in currcentidowns:
			if liner[3] in currboundplus:
				print >> ten, liner[13]
				print >> tenperm, liner[3] + "\t" + liner[13]
				print >> eleven, liner[12]
				print >> elevenperm, liner[3] + "\t" + liner[12]
				print >> twelve, liner[10]
				print >> twelveperm, liner[3] + "\t" + liner[10]
			if liner[3] in currdeplus:
				currcentidowngenes.append(liner[3])
				print >> tenten, liner[13]
				print >> eleveneleven, liner[12]
				print >> twelvetwelve, liner[10]
	master.close()
	if len(currcentifactor) > 0:
		one.close(),oneperm.close(),two.close(),twoperm.close(),three.close(),threeperm.close(),oneone.close(),twotwo.close(),threethree.close()
	if len(currcentidowns) > 0:
		ten.close(),tenperm.close(),eleven.close(),elevenperm.close(),twelve.close(),twelveperm.close(),tenten.close(),eleveneleven.close(),twelvetwelve.close()
	genecounts[gene][2] = len(set(currcentigenes))
	genecounts[gene][3] = len(set(currcentidowngenes))

counting = open(permdir + 'PermCounts.txt','w')
for gener in sorted(genecounts.keys()):
	print >> counting, gener + "\t" + str(genecounts[gener][0]) + "\t" + str(genecounts[gener][1]) + "\t" + str(genecounts[gener][2]) + "\t" + str(genecounts[gener][3])

counting.close()
# print 'Making random estimates...'
# sys.stdout.flush()
# tests = os.listdir(outdir)
# x=1
# for test in sorted(tests):
# #for test in ['GrepChipvCenti']:
# 	if 'Grep' not in test:
# 		continue
# 	masterrandostf = {}
# 	masterrandosdown = {}
# 	for i in range(1000):
# 		masterrandostf[i] = []
# 		masterrandosdown[i] = []
# 	permresults = open(permdir + test + '.master.perms','w')
# 	genes = os.listdir(outdir + test)
# 	for gene in sorted(genes):
# #	for gene in ['YY1_second.down.perm.cage']:
# #		if 'YY1' not in gene:
# #			continue
# 		if '.perm.' not in gene:
# 			continue
# 		if '.tf.' in gene:
# 			currflen = genecounts[gene.split('.')[0]][0]
# 			currclen = genecounts[gene.split('.')[0]][2]
# 			tfscope = 'TF'
# 		if '.down.' in gene:
# 			currflen = genecounts[gene.split('.')[0]][1]
# 			currclen = genecounts[gene.split('.')[0]][3]
# 			tfscope = 'Downstream'
# 		justgene = gene.split('.')[0]
# 		print gene + ' randoms (' + str(currflen) + '/' + str(currclen) + ' genes)...'
# 		sys.stdout.flush()
# 		currdict = {}
# 		currperm = open(outdir + test + '/' + gene,'r')
# 		for line in currperm:
# 			if line == '':
# 				print 'pass'
# 				continue
# 			liner = line.strip().split()
# 			if liner[0] not in currdict.keys():
# 				currdict[liner[0]] = [liner[1]]
# 			else:
# 				currdict[liner[0]].append(liner[1])
# 		currperm.close()
# 		print 'Current dictionary has ' + str(len(currdict)) + ' genes...'
# 		sys.stdout.flush()
# 		#print currdict[currdict.keys()[0]]
# 		generand = []
# 		for g in range(1000):
# 			#print g
# 			random.seed(x)
# 		#	if 'IRF4' in gene:
# 		#		randogenes = sorted(currdict.keys()).extend(random.sample(sorted(currdict.keys()),50))
# 			if 'perm.p' not in gene:
# 				randogenes = random.sample(sorted(currdict.keys()),currflen)
# 			if 'perm.p' in gene:
# 				randogenes = random.sample(sorted(currdict.keys()),currclen)
# 			randos = []
# 			for rando in randogenes:
# 				randos.extend(currdict[rando])
# 			if 'GrepChip' in test:
# 				if '.tf.' in gene:
# 					masterrandostf[g].extend(randos)
# 				if '.down.' in gene:
# 					masterrandosdown[g].extend(randos)
# 				randy = randos.count('centi')/float(len(randos))
# 				generand.append(str(randy))
# 				x += 1
# 				continue
# 			if 'GrepDistance' in test:
# 				randy = [abs(int(z)) for z in randos]
# 				if '.tf.' in gene:
# 					masterrandostf[g].extend(randy)
# 				if '.down.' in gene:
# 					masterrandosdown[g].extend(randy)
# 				generand.append(str(numpy.median(randy)))
# 				x += 1
# 				continue
# 			randy = [float(z) for z in randos]
# 			if '.tf.' in gene:
# 				masterrandostf[g].extend(randy)
# 			if '.down.' in gene:
# 				masterrandosdown[g].extend(randy)
# 			generand.append(str(numpy.median(randy)))
# 			x += 1
# 		print >> permresults, justgene + "\t" + tfscope + "\t" + "\t".join(generand)
# 	permresults.close()
# 	currpermresultstf = open(outdir + test + '.alt.tf.perms','w')
# 	currpermresultsdown = open(outdir + test + '.alt.down.perms','w')
# 	for rep in range(1000):
# 		if 'GrepChip' in test:
# 			print >> currpermresultstf, str(masterrandostf[rep].count('centi')/float(len(masterrandostf[rep])))
# 			print >> currpermresultsdown, str(masterrandosdown[rep].count('centi')/float(len(masterrandosdown[rep])))
# 		else:
# 			print >> currpermresultstf, str(numpy.median(masterrandostf[rep]))
# 			print >> currpermresultsdown, str(numpy.median(masterrandosdown[rep]))
# 	currpermresultstf.close()
# 	currpermresultsdown.close()