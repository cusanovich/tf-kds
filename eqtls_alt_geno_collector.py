#!/usr/bin/env python
import glob

windowsize = 100000
windowname = str(windowsize/1000) + 'kb'
genodir = '/mnt/lustre/home/jdegner/IMPUTED_YRI_7_JUL_11/FINAL.BIMBAM.INPUT/OUTDIR/output/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/'
annots = glob.glob('/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Bindings/*.txt')
exprs = outdir + 'Expression/gene_expression_qqnorm.txt'
indivs = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/gene_expression_levels_noqq.txt').readline().strip().split()
genoindivs = ['pos','minor','major']
genoindivs.extend(open(genodir + 'individual.order.txt','r').readline().strip().strip('IND').split(','))
genocol = [0,1,2]
for indiv in indivs:
	if indiv[0:2] != "NA":
		continue
	genocol.append(genoindivs.index(indiv))

factors = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
factoring = {}
for line in factors:
	if '#' in line:
		continue
	liner = line.strip().split()
	if liner[0] not in factoring.keys():
		factoring[liner[0]] = liner[2]

factors.close()
tss = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/hg18_HT12ensemblTSScombinedsorted.bed','r')
tsss = {}
for line in tss:
	liner = line.strip().split()
	tsss[liner[3]] = liner[0:3]

tss.close()

for annot in annots:
	currannot = annot.split('/')[-1].split('.')[0]
	currclass = annot.split('/')[-1].split('.')[1]
	print currannot + " " + currclass
	exprgenes = open(annot).readline().strip().split('|')
	currens = factoring[currannot]
	currchr = tsss[currens][0]
	if currchr == 'chrX':
		print currannot + ' is on the wrong chromosome!!!'
		continue
	currpos = int(tsss[currens][1])
	print 'Writing expression data...'
	expring = open(exprs,'r')
	curroutexprs = open(outdir + '/AltExpression/' + currclass + '/' + currannot + '_exprs.txt','w')
	curroutexprprime = open(outdir + '/AltExpression/' + currclass + '/' + currannot + '_exprprime.txt','w')
	curroutgenes = open(outdir + '/AltExpression/' + currclass + '/' + currannot + '_geneids.txt','w')
	geneids = expring.readline().strip().split()
	geneind = []
	masterind = geneids.index(currens)
	for gene in exprgenes:
		if gene == currens:
			continue
		if gene in geneids:
			geneind.append(geneids.index(gene))
			print >> curroutgenes, gene
	for line in expring:
		liner = line.strip().split()
		goodexprs = []
		for inder in geneind:
			goodexprs.append(liner[inder])
		print >> curroutexprs, '\t'.join(goodexprs)
		print >> curroutexprprime, liner[masterind]
	curroutexprs.close()
	curroutgenes.close()
	curroutexprprime.close()
	# print 'Writing genotypes...'
	# currinrs = open(genodir + currchr + '.YRI.snpdata.txt','r')
	# rsdic = {}
	# for line in currinrs:
	# 	if 'rs' in line:
	# 		rsdic[line.strip().split()[0]] = line.strip().split()[5]
	# currinrs.close()
	# curringenos = open(genodir + currchr + '.YRI.mean.genotype.txt','r')
	# curroutgenos = open(outdir + '/Genotypes/' + currannot + '_' + windowname + '.mean.genotype.txt','w')
	# for line in curringenos:
	# 	if line == '':
	# 		continue
	# 	if 'rs' in line:
	# 		currdist = int(rsdic[line.strip().split()[0]]) - currpos
	# 	else:
	# 		currdist = int(line.strip().split()[0].split('.')[-1]) - currpos
	# 	if currdist > -windowsize and currdist < windowsize:
	# 		liner = line.strip().split()
	# 		currcalls = []
	# 		for column in genocol:
	# 			currcalls.append(liner[column])
	# 		print >> curroutgenos, '\t'.join(currcalls)
	# curroutgenos.close()