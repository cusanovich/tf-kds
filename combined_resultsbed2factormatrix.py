#!/usr/bin/env python
import sys
from kdfunc import pwmdicter

windowname = '10kb'
resultsbed = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_new_sorted.bed','r')
pwms = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
pwmlist = pwmdicter(pwms)
pwms.close()
genes = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt','r')
resultsmatrix = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/allunionbindingresults' + windowname + '.txt','w')
nopromresultsmatrix = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/noprom_allunionbindingresults' + windowname + '.txt','w')

genelist = []
for line  in genes:
    if "probeID" in line:
        continue
    genelist.append(line.strip().split()[6])

genelist = set(genelist)
genes.close()

factorcount = len(pwmlist)
print factorcount
factororder = sorted(pwmlist.iterkeys())

bindcounts = {}
noprombindcounts = {}
for gene in genelist:
    bindcounts[gene] = [0]*factorcount
    noprombindcounts[gene] = [0]*factorcount

print 'Counting binding events...'
#bindcounts = {}
for line in resultsbed:
    geneid = line.strip().split()[3]
    if len(line.strip().split()) > 9 and geneid in genelist:
        binder = line.strip().split()[9]
        if binder[0:2] == "M0" or "_" in binder:
            for factor in factororder:
                if binder in pwmlist[factor]:
                    bindcounts[geneid][factororder.index(factor)] += 1
                    if abs(int(line.strip().split()[1]) - int(line.strip().split()[7])) > 1000:
                        noprombindcounts[geneid][factororder.index(factor)] += 1
        else:
            bindcounts[geneid][factororder.index(binder)] += 1
            if abs(int(line.strip().split()[1]) - int(line.strip().split()[7])) > 1000:
                noprombindcounts[geneid][factororder.index(binder)] += 1

print 'Writing results...'
resultsbed.close()
print >> resultsmatrix, '"' + '"\t"'.join(factororder) + '"'
print >> nopromresultsmatrix, '"' + '"\t"'.join(factororder) + '"'
for gene in sorted(bindcounts.iterkeys()):
    genic = [str(x) for x in bindcounts[gene]]
    print >> resultsmatrix, '"' + gene + '"\t' + '\t'.join(genic)
    nopromgenic = [str(x) for x in noprombindcounts[gene]]
    print >> nopromresultsmatrix, '"' + gene + '"\t' + '\t'.join(nopromgenic)

resultsmatrix.close()
nopromresultsmatrix.close()