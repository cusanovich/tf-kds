#!/usr/bin/env python
import sys

windowname = '10kb'
resultsbed = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_sorted.bed','r')
pwms = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
resultsmatrix = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/allbindingresults' + windowname + '.txt','w')

pwmlist = {}
for line in pwms:
    if '#' in line:
        continue
    if line.strip().split()[0] in pwmlist.keys():
        pwmlist[line.strip().split()[0]].append(line.strip().split()[1])
    if line.strip().split()[0] not in pwmlist.keys():
        pwmlist[line.strip().split()[0]] = [line.strip().split()[1]]

pwms.close()
factorcount = len(pwmlist)
print factorcount
factororder = sorted(pwmlist.iterkeys())

print 'Counting binding events...'
print >> resultsmatrix, '"' + '"\t"'.join(factororder) + '"'
bindcounts = {}
for line in resultsbed:
    geneid = line.strip().split()[3]
    if geneid not in bindcounts.keys():
        bindcounts[geneid] = [0]*factorcount
    if len(line.strip().split()) > 9:
        binder = line.strip().split()[9]
        for factor in factororder:
            if binder in pwmlist[factor]:
                bindcounts[geneid][factororder.index(factor)] += 1

print 'Writing results...'
resultsbed.close()
for gene in sorted(bindcounts.iterkeys()):
    genic = [str(x) for x in bindcounts[gene]]
    print >> resultsmatrix, '"' + gene + '"\t' + '\t'.join(genic)

resultsmatrix.close()