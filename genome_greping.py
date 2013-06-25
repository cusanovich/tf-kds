#!/usr/bin/env python
import os

grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/'
genelists = os.listdir(grepdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))
#genelist = ['ARNTL2_old']
factors = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt'
masterfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/results10kb_combined_sorted.bed'
centifile = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/10kbresults_phastcons.bed'

def genereader(thisfile):
	currfile = open(thisfile,'r')
	return(currfile.readline().strip().split('|'))

for gene in genelist:
	print gene
	currtf = genereader(grepdir+gene+'.BindingTF.txt')
	currdownstream = genereader(grepdir+gene+'.DownstreamTFs.txt')
	currbound = genereader(grepdir+gene+'.Bound.txt')
	currboundplus = genereader(grepdir+gene+'.BoundPlus.txt')
	currde = genereader(grepdir+gene+'.DEandBound.txt')
	currdeplus = genereader(grepdir+gene+'.DEandBoundPlus.txt')
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

	print 'Splitting Master Records...'
	master = open(masterfile,'r')
	tfcagegenes = []
	downcagegenes = []
	if len(currfactor) > 0:
		one = open(outdir + 'GrepDistance/' + gene + '.tf.bound.distance','w')
		two = open(outdir + 'GrepCage/' + gene + '.tf.bound.cage','w')
		three = open(outdir + 'GrepChipvCenti/' + gene + '.tf.bound.chipvcenti','w')
		oneone = open(outdir + 'GrepDistance/' + gene + '.tf.debound.distance','w')
		twotwo = open(outdir + 'GrepCage/' + gene + '.tf.debound.cage','w')
		threethree = open(outdir + 'GrepChipvCenti/' + gene + '.tf.debound.chipvcenti','w')
	if len(currdowns) > 0:
		ten = open(outdir + 'GrepDistance/' + gene + '.down.boundplus.distance','w')
		eleven = open(outdir + 'GrepCage/' + gene + '.down.boundplus.cage','w')
		twelve = open(outdir + 'GrepChipvCenti/' + gene + '.down.boundplus.chipvcenti','w')
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
				print >> one, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))*strand
				if liner[3] not in tfcagegenes:
					if '_' in liner[4]:
						print >> two, liner[4].split('_')[0]
					else:
						print >> two, liner[4]
					tfcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> three, chiper
			if liner[3] in currde:
				print >> oneone, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))*strand
				if liner[3] not in tfcagegenes:
					if '_' in liner[4]:
						print >> twotwo, liner[4].split('_')[0]
					else:
						print >> twotwo, liner[4]
					tfcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> threethree, chiper
		if liner[9] in currdowns:
			if liner[3] in currboundplus:
				print >> ten, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))*strand
				if liner[3] not in downcagegenes:
					if '_' in liner[4]:
						print >> eleven, liner[4].split('_')[0]
					else:
						print >> eleven, liner[4]
					downcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> twelve, chiper
			if liner[3] in currdeplus:
				print >> tenten, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))*strand
				if liner[3] not in downcagegenes:
					if '_' in liner[4]:
						print >> eleveneleven, liner[4].split('_')[0]
					else:
						print >> eleveneleven, liner[4]
					downcagegenes.append(liner[3])
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> twelvetwelve, chiper

	if len(currfactor) > 0:
		one.close(),two.close(),three.close(),oneone.close(),twotwo.close(),threethree.close()
	if len(currdowns) > 0:
		ten.close(),eleven.close(),twelve.close(),tenten.close(),eleveneleven.close(),twelvetwelve.close()

	print 'Splitting Centipede Records...'
	master = open(centifile,'r')
	if len(currcentifactor) > 0:
		one = open(outdir + 'GrepPhastCons/' + gene + '.tf.bound.phast','w')
		two = open(outdir + 'GrepPWM/' + gene + '.tf.bound.pwm','w')
		three = open(outdir + 'GrepPosterior/' + gene + '.tf.bound.post','w')
		oneone = open(outdir + 'GrepPhastCons/' + gene + '.tf.debound.phast','w')
		twotwo = open(outdir + 'GrepPWM/' + gene + '.tf.debound.pwm','w')
		threethree = open(outdir + 'GrepPosterior/' + gene + '.tf.debound.post','w')
	if len(currcentidowns) > 0:
		ten = open(outdir + 'GrepPhastCons/' + gene + '.down.boundplus.phast','w')
		eleven = open(outdir + 'GrepPWM/' + gene + '.down.boundplus.pwm','w')
		twelve = open(outdir + 'GrepPosterior/' + gene + '.down.boundplus.post','w')
		tenten = open(outdir + 'GrepPhastCons/' + gene + '.down.deboundplus.phast','w')
		eleveneleven = open(outdir + 'GrepPWM/' + gene + '.down.deboundplus.pwm','w')
		twelvetwelve = open(outdir + 'GrepPosterior/' + gene + '.down.deboundplus.post','w')
	for line in master:
		liner = line.strip().split()
		if liner[9] in currcentifactor:
			if liner[3] in currbound:
				print >> one, liner[13]
				print >> two, liner[12]
				print >> three, liner[10]
			if liner[3] in currde:
				print >> oneone, liner[13]
				print >> twotwo, liner[12]
				print >> threethree, liner[10]
		if liner[9] in currcentidowns:
			if liner[3] in currboundplus:
				print >> ten, liner[13]
				print >> eleven, liner[12]
				print >> twelve, liner[10]
			if liner[3] in currdeplus:
				print >> tenten, liner[13]
				print >> eleveneleven, liner[12]
				print >> twelvetwelve, liner[10]

	if len(currcentifactor) > 0:
		one.close(),two.close(),three.close(),oneone.close(),twotwo.close(),threethree.close()
	if len(currcentidowns) > 0:
		ten.close(),eleven.close(),twelve.close(),tenten.close(),eleveneleven.close(),twelvetwelve.close()
