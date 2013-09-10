#!/usr/bin/env python
import time
import subprocess
import numpy
import random
import sys
import linecache

uniqid = str(int(time.time()*100000))
perms = 9000
gene = sys.argv[1]
#gene = 'IRF3'
masterexprdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/'
exprdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/Bound/'
genodir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Genotypes/'
genos = genodir + gene + '_100kb.mean.genotype.txt'
primexpr = exprdir + gene + '_exprprime.txt'
phenos = exprdir + gene + '_exprs.txt'
mastercovs = masterexprdir + 'gene_expression_cov.txt'
geneids = open(exprdir + gene + '_geneids.txt').readlines()
geneids = [x.strip().split()[0] for x in geneids]
express = numpy.loadtxt(exprdir + gene + '_exprs.txt')
primexpress = numpy.loadtxt(primexpr)
exprs = express.shape[1]
cover = numpy.loadtxt(masterexprdir + 'gene_expression_cov.txt')
covsize = cover.shape[0]
resultant = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/eQTLs/Bound/' + gene + '_100kb_sidney_results.txt'
reupdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/eQTLs/Bound/Reups/'
permexpr = './output/' + uniqid + '.temp'
permcov = './output/' + uniqid + '.cov'

def gemming(genos=0,phenos=0,covs=0,uniqid=0,i=0):
	gemmer = '/mnt/lustre/home/cusanovich/Programs/gemma -g ' + genos + ' -p ' + phenos + ' -k ' + covs + ' -lmm 1 -n ' + str(i) + ' -o ' + uniqid
	gemmify = subprocess.Popen(gemmer,shell=True)
	gemmify.wait()
	dt = numpy.dtype('S1,S50,S1,S1,S1,S1,S1,f4')
	report = numpy.loadtxt('./output/' + uniqid + '.assoc.txt',delimiter='\t',dtype=dt,skiprows=1)
	try:
		pvals = [x[7] for x in report]
		pval = min(pvals)
		rss = [x[1] for x in report]
		rs = rss[pvals.index(pval)]
	except TypeError:
		pval = float(linecache.getline('./output/' + uniqid + '.assoc.txt', 2).strip().split()[7])
		rs = linecache.getline('./output/' + uniqid + '.assoc.txt', 2).strip().split()[1]	
		linecache.clearcache()
	return([rs,pval])

checker = open(resultant,'r')
reups = {}
x=0
for line in checker:
	liner = line.strip().split()
	if 'self-qtl' in line:
		masterqtl = liner
	else:
		reups[x] = liner
		x += 1
checker.close()

currindex = range(covsize)
if float(masterqtl[4])*1000 < 10:
	mastercalib = int(float(masterqtl[4])*1000)
	for m in range(perms):
		random.shuffle(currindex)
		currexpress = primexpress[currindex,:]
		currcov = cover[currindex,:]
		currcov = currcov[:,currindex]
		numpy.savetxt(permexpr,currexpress,delimiter='\t')
		numpy.savetxt(permcov,currcov,delimiter='\t')
		permqtl = gemming(genos=genos,phenos=permexpr,covs=permcov,uniqid=uniqid,i=1)
		if permqtl[1] < float(masterqtl[3]):
			mastercalib += 1
	masterqtl[4] = str(float(mastercalib)/10000.0)

reupper = open(reupdir + gene + '_100kb_sidney_reups.txt','w')
print >> reupper, '\t'.join(masterqtl)
genfilter = open(genos,'r')
genoholder = './output/' + uniqid + '.mean.genotype.txt'
genoqtl = open(genoholder,'w')
for line in genfilter:
	if masterqtl[1] in line:
		print >> genoqtl, line.strip()

genoqtl.close()

for k in sorted(reups.keys()):
#for k in [0]:
	currreups = reups[k]
	if float(currreups[4])*1000 < 10:
		permp = int(float(currreups[4])*1000)
		for j in range(perms):
			random.shuffle(currindex)
			currexpress = express[currindex,:]
			currcov = cover[currindex,:]
			currcov = currcov[:,currindex]
			numpy.savetxt(permexpr,currexpress[:,k],delimiter='\t')
			numpy.savetxt(permcov,currcov,delimiter='\t')
			permer = gemming(genos=genoholder,phenos=permexpr,covs=permcov,uniqid=uniqid,i=1)
			if permer[1] < float(currreups[3]):
				permp += 1
		currreups[4] = str(float(permp)/10000.0)
	print >> reupper, '\t'.join(currreups)

reupper.close()