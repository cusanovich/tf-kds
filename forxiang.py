#!/usr/bin/env python2.6
import os
import sys
import subprocess
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
from pybedtools import BedTool, featurefuncs
from kdfunc import genereader
#I want a dictionary with each motif linked to functionally bound genes
#I need to know if a motif is binding a gene and whether that gene is DE in any of the appropriate knockdowns

windowsize = '10kb'
grepdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/Union/' + windowsize + '/'
factors = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
centifile = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'results_phastcons_sorted.bed','r')
targetfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'results_sorted_forxiang.bed'
targetplusfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'resultsplus_sorted_forxiang.bed'

#This establishes a dictionary where each motif is linked to a gene symbol (so that the target dictionaries can be converted)
factordic = {}
for line in factors:
	liner = line.strip().split()
	if '_' in liner[1]:
		continue
	if liner[0] not in factordic.keys():
		factordic[liner[0]] = [liner[1]]
	else:
		factordic[liner[0]].append(liner[1])

factors.close()

genelists = os.listdir(grepdir)
genelist = sorted(set([x.split('.')[0] for x in genelists]))

targets = {}
targetsplus = {}
for gene in genelist:
	print gene
	if gene+'.DEandBound.txt' in genelists:
		currtargets = genereader(grepdir+gene+'.DEandBound.txt')
		currfactor = genereader(grepdir+gene+'.BindingTF.txt')[0]
		try:
			currmotifs = factordic[currfactor]
		except KeyError:
			pass
		for motif in currmotifs:
			if motif not in targets.keys():
				targets[motif] = currtargets
			if motif in targets.keys():
				temp = targets[motif]
				targets[motif] = list(set(temp+currtargets))
	currtargetsplus = genereader(grepdir+gene+'.DEandBoundPlus.txt')
	currfactors = genereader(grepdir+gene+'.DownstreamTFs.txt')
	for factor in currfactors:
		try:
			currmotifs = factordic[factor]
		except KeyError:
			continue
		for motif in currmotifs:
			if motif not in targetsplus.keys():
				targetsplus[motif] = currtargetsplus
			if motif in targetsplus.keys():
				temp = targetsplus[motif]
				targetsplus[motif] = list(set(temp+currtargetsplus))

print 'Collecting functional associations...'
funcbinds = {}
funcbindsplus = {}
for line in centifile:
	liner = line.strip().split()
	if len(liner) < 13:
		continue
	adder = ["\t".join(liner[6:9]) + "\t" + liner[3]]
	if liner[9] in targets and liner[3] in targets[liner[9]]:
		try:
			funcbinds[liner[9]] = funcbinds[liner[9]] + adder
		except KeyError:
			funcbinds[liner[9]] = adder
	if liner[9] in targetsplus and liner[3] in targetsplus[liner[9]]:
		try:
			adder = ["\t".join(liner[6:9]) + "\t" + liner[3]]
			funcbindsplus[liner[9]] = funcbindsplus[liner[9]] + adder
		except KeyError:
			funcbindsplus[liner[9]] = adder

centifile.close()

def name2score(f):
    f = featurefuncs.extend_fields(f, 5)
    if ';' in f.name:
    	f.name = ";".join(list(sorted(set(f.name.split(';')))))
    f.score, f.name = f.name, motifer
    return f

print 'Printing direct targets...'
for motifer in funcbinds.keys():
	print motifer
	currmotif = BedTool(funcbinds[motifer]).sort().merge(nms=True).each(name2score).saveas('./temp1')
	cleaner = 'cat ./temp1 >> ./temp2; rm ./temp1'
	cleaning = subprocess.Popen(cleaner,shell=True)
	cleaning.wait()

sorter = 'grep . ./temp2 | sort -k1,1 -k2,2n -k3,3n -k4,4 > ' + targetfile + '; rm ./temp2'
sorting = subprocess.Popen(sorter,shell=True)
sorting.wait

print 'Printing indirect targets...'
for motifer in funcbindsplus.keys():
	print motifer
	currmotif = BedTool(funcbindsplus[motifer]).sort().merge(nms=True).each(name2score).saveas('./temp1')
	cleaner = 'cat ./temp1 >> ./temp2; rm ./temp1'
	cleaning = subprocess.Popen(cleaner,shell=True)
	cleaning.wait()

sorter = 'grep . ./temp2 | sort -k1,1 -k2,2n -k3,3n -k4,4 > ' + targetplusfile + '; rm ./temp2'
sorting = subprocess.Popen(sorter,shell=True)
sorting.wait