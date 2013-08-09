#!/usr/bin/env python

windowname = '10kb'
pwms = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
states = open('/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/chromstates_list.txt','r')
pather = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Matrices/' + windowname + '/'
genes = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt','r')
olap = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/' + windowname + '_allbinding_chromstatesandtss_overlap_sorted.bed'

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

statelist = states.readlines()
statelist = [x.strip() for x in statelist]
states.close()

genelist = genes.readlines()
genelist = [x.strip().split()[6] for x in genelist]
genes.close()

for state in statelist:
    print state
    resultsbed = open(olap,'r')
    resultsmatrix = open(pather + state + "_resultsmatrix.txt",'w')
    print 'Counting binding events...'
    bindcounts = {}
    for line in resultsbed:
        if len(line.strip().split()) > 9 and state not in line.strip().split()[14]:
            continue
        if line.strip().split()[0] == 'chrY':
            continue
        geneid = line.strip().split()[3]
        if geneid not in bindcounts.keys():
            bindcounts[geneid] = [0]*factorcount
        if len(line.strip().split()) > 9 and state in line.strip().split()[14]:
            binder = line.strip().split()[9]
            for factor in factororder:
                if binder in pwmlist[factor]:
                    bindcounts[geneid][factororder.index(factor)] += 1

    print 'Writing results...'
    print >> resultsmatrix, '"' + '"\t"'.join(factororder) + '"'
    placeholder = [0]*factorcount
    placeholder = [str(x) for x in placeholder]
    #for gene in sorted(bindcounts.iterkeys()):
    for gene in sorted(genelist):
        if gene in bindcounts.iterkeys():
            genic = [str(x) for x in bindcounts[gene]]
            print >> resultsmatrix, '"' + gene + '"\t' + '\t'.join(genic)
        else:
            print >> resultsmatrix, '"' + gene + '"\t' + '\t'.join(placeholder)

    resultsmatrix.close()
    resultsbed.close()