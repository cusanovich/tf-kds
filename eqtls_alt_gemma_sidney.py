#!/usr/bin/env python
import time
import subprocess
import numpy
import random
import sys
import linecache

uniqid = str(int(time.time()*100000))
perms = 1000
gene = sys.argv[1]
#gene = 'IRF3'
#bindingclass = 'DEandBound'
bindingclass = sys.argv[2]
masterexprdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/'
exprdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltExpression/' + bindingclass + '/'
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
currindex = range(covsize)
resultant = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AlteQTLs/' + bindingclass + '/' + gene + '_100kb_sidney_results.txt','w')

permexpr = './output/' + uniqid + '.temp'
permcov = './output/' + uniqid + '.cov'

def gemming(genos=0,phenos=0,covs=0,uniqid=uniqid,i=1):
	gemmer = '/mnt/lustre/home/cusanovich/Programs/gemma -g ' + genos + ' -p ' + phenos + ' -k ' + covs + ' -lmm 1 -n ' + str(i) + ' -maf 0.05 -o ' + str(uniqid)
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

def perming(obser,express,genos,permp=0,perms=perms,index=currindex,cover=cover,permexpr=permexpr,permcov=permcov,k=0):
	for j in range(perms):
		print "..........\n" + str(j) + "\n........."
		random.shuffle(currindex)
		currexpress = express[currindex,:]
		currcov = cover[currindex,:]
		currcov = currcov[:,currindex]
		try:
			currexpress.shape[1]
			numpy.savetxt(permexpr,currexpress[:,k],delimiter='\t')
		except IndexError:
			numpy.savetxt(permexpr,currexpress,delimiter='\t')
		numpy.savetxt(permcov,currcov,delimiter='\t')
		permer = gemming(genos=genos,phenos=permexpr,covs=permcov)
		if permer[1] < obser[1]:
			permp += 1
	return(permp)

masterqtl = gemming(genos=genos,phenos=primexpr,covs=mastercovs)
runperms = perms
mastercalib = perming(obser=masterqtl,express=primexpress,genos=genos)
if mastercalib >= 50:
	masterqtl.append((float(mastercalib)/float(runperms)))
	winner = [gene,masterqtl[0],'self-qtl',str(masterqtl[1]),str(masterqtl[2])]
	print >> resultant, '\t'.join(winner)
	sys.exit(0)
if mastercalib < 50 and bindingclass == "DEandBound":
	masterqtl.append((float(mastercalib)/float(runperms)))
	winner = [gene,masterqtl[0],'self-qtl',str(masterqtl[1]),str(masterqtl[2])]
	reportant = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/RealCis/' + gene + '_100kb_sidney_results.txt','w')
	print >> reportant, '\t'.join(winner)
	reportant.close()
if mastercalib < 10:
	mastercalib = perming(obser=masterqtl,express=primexpress,genos=genos,permp=mastercalib,perms=9000)
	runperms = 10000

masterqtl.append((float(mastercalib)/float(runperms)))
winner = [gene,masterqtl[0],'self-qtl',str(masterqtl[1]),str(masterqtl[2])]
print >> resultant, '\t'.join(winner)
genfilter = open(genos,'r')
genoholder = './output/' + uniqid + '.mean.genotype.txt'
genoqtl = open(genoholder,'w')
for line in genfilter:
	if masterqtl[0] in line:
		print >> genoqtl, line.strip()

genoqtl.close()

for k in range(exprs):
	print "---------\n" + str(k) + " of " + str(exprs) + "\n---------"
#for k in [0]:
	runperms = perms
	print 'Running GEMMA...'
	transer = gemming(genos=genoholder,phenos=phenos,covs=mastercovs,i=(k+1))
	print 'Starting permutations...'
	permp = perming(obser=transer,express=express,genos=genoholder)
	if permp < 10:
		permp = perming(obser=transer,express=express,genos=genoholder,permp=permp,perms=9000)
		runperms = 10000

	winner = [gene,masterqtl[0],geneids[k],str(transer[1]),str(float(permp)/float(runperms))]
	print >> resultant, '\t'.join(winner)

resultant.close()