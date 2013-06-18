#!/usr/bin/env python
import os

#if(len(sys.argv)!= 4):
#	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: gff2bed.py [Outdir] [In File] [Out File]\n\n')
#	sys.exit(1)

#motif = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/all_motifs.txt','r')
#motifs = motif.readline().strip().split('|')
grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/'
genelists = os.listdir(grepdir)
#phenolists = os.listdir(outdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))
#genelist = ['ARNTL2_old']
#phenolist = sorted(set([x.split('Grep')[1] for x in phenolists]))
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
	for line in currf:
		if line.strip().split()[0] in currtf:
			currfactor.append(line.strip().split()[1])
		if line.strip().split()[0] in currdownstream:
			currdowns.append(line.strip().split()[1])

	currf.close()
	currfactor = set(currfactor)
	currdowns = set(currdowns)

	print 'Splitting Master Records...'
	master = open(masterfile,'r')
	if len(currfactor) > 0:
		one = open(outdir + 'GrepDistance/' + gene + '.tf.bound.distance','w')
		two = open(outdir + 'GrepCage/' + gene + '.tf.bound.cage','w')
		three = open(outdir + 'GrepChipvCenti/' + gene + '.tf.bound.chipvcenti','w')
		four = open(outdir + 'GrepDistance/' + gene + '.tf.boundplus.distance','w')
		five = open(outdir + 'GrepCage/' + gene + '.tf.boundplus.cage','w')
		six = open(outdir + 'GrepChipvCenti/' + gene + '.tf.boundplus.chipvcenti','w')
		oneone = open(outdir + 'GrepDistance/' + gene + '.tf.debound.distance','w')
		twotwo = open(outdir + 'GrepCage/' + gene + '.tf.debound.cage','w')
		threethree = open(outdir + 'GrepChipvCenti/' + gene + '.tf.debound.chipvcenti','w')
		fourfour = open(outdir + 'GrepDistance/' + gene + '.tf.deboundplus.distance','w')
		fivefive = open(outdir + 'GrepCage/' + gene + '.tf.deboundplus.cage','w')
		sixsix = open(outdir + 'GrepChipvCenti/' + gene + '.tf.deboundplus.chipvcenti','w')
	if len(currdowns) > 0:
		seven = open(outdir + 'GrepDistance/' + gene + '.down.bound.distance','w')
		eight = open(outdir + 'GrepCage/' + gene + '.down.bound.cage','w')
		nine = open(outdir + 'GrepChipvCenti/' + gene + '.down.bound.chipvcenti','w')
		ten = open(outdir + 'GrepDistance/' + gene + '.down.boundplus.distance','w')
		eleven = open(outdir + 'GrepCage/' + gene + '.down.boundplus.cage','w')
		twelve = open(outdir + 'GrepChipvCenti/' + gene + '.down.boundplus.chipvcenti','w')
		sevenseven = open(outdir + 'GrepDistance/' + gene + '.down.debound.distance','w')
		eighteight = open(outdir + 'GrepCage/' + gene + '.down.debound.cage','w')
		ninenine = open(outdir + 'GrepChipvCenti/' + gene + '.down.debound.chipvcenti','w')
		tenten = open(outdir + 'GrepDistance/' + gene + '.down.deboundplus.distance','w')
		eleveneleven = open(outdir + 'GrepCage/' + gene + '.down.deboundplus.cage','w')
		twelvetwelve = open(outdir + 'GrepChipvCenti/' + gene + '.down.deboundplus.chipvcenti','w')
	for line in master:
		liner = line.strip().split()
		if len(liner) < 7:
			continue
		if liner[9] in currfactor:
			if liner[3] in currbound:
				print >> one, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> two, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> three, chiper
			if liner[3] in currboundplus:
				print >> four, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> five, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> six, chiper
				continue
			if liner[3] in currde:
				print >> oneone, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> twotwo, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> threethree, chiper
			if liner[3] in currdeplus:
				print >> fourfour, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> fivefive, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> sixsix, chiper
				continue
		if liner[9] in currdowns:
			if liner[3] in currbound:
				print >> seven, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> eight, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> nine, chiper
			if liner[3] in currboundplus:
				print >> ten, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> eleven, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> twelve, chiper
				continue
			if liner[3] in currde:
				print >> sevenseven, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> eighteight, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> ninenine, chiper
			if liner[3] in currdeplus:
				print >> tenten, int(round(((float(liner[7])+float(liner[8]))/2) - float(liner[2]),0))
				print >> eleveneleven, liner[4]
				chiper = 'chip'
				if len(liner[9].split('_')) < 2:
					chiper = 'centi'
				print >> twelvetwelve, chiper
				continue

	if len(currfactor) > 0:
		one.close(),two.close(),three.close(),four.close(),five.close(),six.close(),oneone.close(),twotwo.close(),threethree.close(),fourfour.close(),fivefive.close(),sixsix.close()
	if len(currdowns) > 0:
		seven.close(),eight.close(),nine.close(),ten.close(),eleven.close(),twelve.close(),sevenseven.close(),eighteight.close(),ninenine.close(),tenten.close(),eleveneleven.close(),twelvetwelve.close()

	print 'Splitting Centipede Records...'
	master = open(centifile,'r')
	if len(currfactor) > 0:
		one = open(outdir + 'GrepPhastCons/' + gene + '.tf.bound.phast','w')
		two = open(outdir + 'GrepPWM/' + gene + '.tf.bound.pwm','w')
		three = open(outdir + 'GrepPosterior/' + gene + '.tf.bound.post','w')
		four = open(outdir + 'GrepPhastCons/' + gene + '.tf.boundplus.phast','w')
		five = open(outdir + 'GrepPWM/' + gene + '.tf.boundplus.pwm','w')
		six = open(outdir + 'GrepPosterior/' + gene + '.tf.boundplus.post','w')
		oneone = open(outdir + 'GrepPhastCons/' + gene + '.tf.debound.phast','w')
		twotwo = open(outdir + 'GrepPWM/' + gene + '.tf.debound.pwm','w')
		threethree = open(outdir + 'GrepPosterior/' + gene + '.tf.debound.post','w')
		fourfour = open(outdir + 'GrepPhastCons/' + gene + '.tf.deboundplus.phast','w')
		fivefive = open(outdir + 'GrepPWM/' + gene + '.tf.deboundplus.pwm','w')
		sixsix = open(outdir + 'GrepPosterior/' + gene + '.tf.deboundplus.post','w')
	if len(currdowns) > 0:
		seven = open(outdir + 'GrepPhastCons/' + gene + '.down.bound.phast','w')
		eight = open(outdir + 'GrepPWM/' + gene + '.down.bound.pwm','w')
		nine = open(outdir + 'GrepPosterior/' + gene + '.down.bound.post','w')
		ten = open(outdir + 'GrepPhastCons/' + gene + '.down.boundplus.phast','w')
		eleven = open(outdir + 'GrepPWM/' + gene + '.down.boundplus.pwm','w')
		twelve = open(outdir + 'GrepPosterior/' + gene + '.down.boundplus.post','w')
		sevenseven = open(outdir + 'GrepPhastCons/' + gene + '.down.debound.phast','w')
		eighteight = open(outdir + 'GrepPWM/' + gene + '.down.debound.pwm','w')
		ninenine = open(outdir + 'GrepPosterior/' + gene + '.down.debound.post','w')
		tenten = open(outdir + 'GrepPhastCons/' + gene + '.down.deboundplus.phast','w')
		eleveneleven = open(outdir + 'GrepPWM/' + gene + '.down.deboundplus.pwm','w')
		twelvetwelve = open(outdir + 'GrepPosterior/' + gene + '.down.deboundplus.post','w')
	for line in master:
		liner = line.strip().split()
		if liner[9] in currfactor:
			if liner[3] in currbound:
				print >> one, liner[13]
				print >> two, liner[12]
				print >> three, liner[10]
			if liner[3] in currboundplus:
				print >> four, liner[13]
				print >> five, liner[12]
				print >> six, liner[10]
				continue
			if liner[3] in currde:
				print >> oneone, liner[13]
				print >> twotwo, liner[12]
				print >> threethree, liner[10]
			if liner[3] in currdeplus:
				print >> fourfour, liner[13]
				print >> fivefive, liner[12]
				print >> sixsix, liner[10]
				continue
		if liner[9] in currdowns:
			if liner[3] in currbound:
				print >> seven, liner[13]
				print >> eight, liner[12]
				print >> nine, liner[10]
			if liner[3] in currboundplus:
				print >> ten, liner[13]
				print >> eleven, liner[12]
				print >> twelve, liner[10]
				continue
			if liner[3] in currde:
				print >> sevenseven, liner[13]
				print >> eighteight, liner[12]
				print >> ninenine, liner[10]
			if liner[3] in currdeplus:
				print >> tenten, liner[13]
				print >> eleveneleven, liner[12]
				print >> twelvetwelve, liner[10]
				continue

	if len(currfactor) > 0:
		one.close(),two.close(),three.close(),four.close(),five.close(),six.close(),oneone.close(),twotwo.close(),threethree.close(),fourfour.close(),fivefive.close(),sixsix.close()
	if len(currdowns) > 0:
		seven.close(),eight.close(),nine.close(),ten.close(),eleven.close(),twelve.close(),sevenseven.close(),eighteight.close(),ninenine.close(),tenten.close(),eleveneleven.close(),twelvetwelve.close()
