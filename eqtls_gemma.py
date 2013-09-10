#!/usr/bin/env python
import time
import subprocess
import numpy
import random
import sys

uniqid = str(int(time.time()*100000))
perms = 1000
gene = sys.argv[1]
#gene = 'ARNTL2'
exprdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/'
genodir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Genotypes/'
exprseek = exprdir + gene + '_exprs.txt'
geneseek = exprdir + gene + '_geneids.txt'
geneids = open(geneseek).readlines()
geneids = [x.strip().split()[0] for x in geneids]
express = numpy.loadtxt(exprseek)
exprs = express.shape[1]
cover = numpy.loadtxt(exprdir + 'gene_expression_cov.txt')
covsize = cover.shape[0]
resultant = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/eQTLs/' + gene + '_results.txt','w')

for i in range(exprs):
	print 'Running GEMMA...'
	gemmer = '/mnt/lustre/home/cusanovich/Programs/gemma -g /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Genotypes/' + gene + '_10kb.mean.genotype.txt -p /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/' + gene + '_exprs.txt -k /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/gene_expression_cov.txt -lmm 1 -n ' + str(i+1) + ' -o ' + uniqid
	gemmify = subprocess.Popen(gemmer,shell=True)
	gemmify.wait()
	dt = numpy.dtype('S1,S50,S1,S1,S1,S1,S1,f4')
	report = numpy.loadtxt('./output/' + uniqid + '.assoc.txt',delimiter='\t',dtype=dt,skiprows=1)
	pvals = [x[7] for x in report]
	pval = min(pvals)
	rss = [x[1] for x in report]
	rs = rss[pvals.index(pval)]
	permp = 0
	print 'Starting permutations...'
	for j in range(perms):
		currindex = range(covsize)
		random.shuffle(currindex)
		currexpress = express[currindex,:]
		currcov = cover[currindex,:]
		currcov = currcov[:,currindex]
		numpy.savetxt('./output/' + uniqid + '.temp',currexpress[:,i],delimiter='\t')
		numpy.savetxt('./output/' + uniqid + '.cov',currcov,delimiter='\t')
		gemmer = '/mnt/lustre/home/cusanovich/Programs/gemma -g /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Genotypes/' + gene + '_10kb.mean.genotype.txt -p ./output/' + uniqid + '.temp -k ./output/' + uniqid + '.cov -lmm 1 -o temp_' + uniqid
		gemmify = subprocess.Popen(gemmer,shell=True)
		gemmify.wait()
		dt = numpy.dtype('S1,S50,S1,S1,S1,S1,S1,f4')
		permer = numpy.loadtxt('./output/temp_' + uniqid + '.assoc.txt',delimiter='\t',dtype=dt,skiprows=1)
		permps = [x[7] for x in permer]
		permper = min(permps)
		if permper < pval:
			permp += 1
	winner = [gene,rs,geneids[i],str(pval),str(float(permp)/float(perms))]
	print >> resultant, '\t'.join(winner)
