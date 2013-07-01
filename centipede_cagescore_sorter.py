#!/usr/bin/env python
import subprocess

print 'Converting GTF to BED...'
cager = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/wgEncodeRikenCageGm12878CellPapTssGencV7.gtf','r')
outter = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSS2.bed','w')
for line in cager:
    liner = line.strip().split()
    chrom = liner[0]
    start = str(int(liner[3]) - 1)
    end = liner[4]
    score = liner[5]
    strand = liner[6]
    gene = liner[9].strip('"').split('.')[0]
    biotype = liner[13].strip('"').split(',')[0]
    print >> outter, chrom + '\t' + start + '\t' + end + '\t' + gene + '\t' + score + '\t' + strand + '\t' + biotype

cager.close()
outter.close()

print 'Consolidating genes...'
bed = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSS2.bed','r')
newbedfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombined2.bed'
newbed = open(newbedfile,'w')
sortbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted2.bed'
print 'Building gene dictionary...'
chrs = {}
starts = {}
ends = {}
scores = {}
strands = {}
biotypes = {}

for line in bed:
    if line == '':
        continue
    liner = line.strip().split()
    try:
        chrs[liner[3]].append(liner[0])
    except KeyError:
        chrs[liner[3]] = [liner[0]]
    try:
        starts[liner[3]].append(int(liner[1]))
    except KeyError:
        starts[liner[3]] = [int(liner[1])]
    try:
        ends[liner[3]].append(int(liner[2]))
    except KeyError:
        ends[liner[3]] = [int(liner[2])]
    try:
        scores[liner[3]].append(float(liner[4]))
    except KeyError:
        scores[liner[3]] = [float(liner[4])]
    try:
        strands[liner[3]].append(liner[5])
    except KeyError:
        strands[liner[3]] = [liner[5]]
    try:
        biotypes[liner[3]].append(liner[6])
    except KeyError:
        biotypes[liner[3]] = [liner[6]]

bed.close()

print 'Dictionary contains ' + str(len(chrs.keys())) + ' genes...'
print 'Collapsing records...'
j=0
k=0
m=0
for bedrecord in chrs.keys():
    if len(chrs[bedrecord]) == 1:
        print >> newbed, chrs[bedrecord][0] + '\t' + str(starts[bedrecord][0]) + '\t' + str(ends[bedrecord][0]) + '\t' + bedrecord + '\t' + str(scores[bedrecord][0]) + '\t' + strands[bedrecord][0]
        j += 1
        continue
    chra = chrs[bedrecord][0]
    stranda = strands[bedrecord][0]
    if biotypes[bedrecord].count('protein_coding') == 0:
        winning = max(scores[bedrecord])
        winners = []
        count = -1
        for score in scores[bedrecord]:
            count += 1
            if score == winning:
                winners.append(count)
        newstarts = []
        wins = [str(winning)]
        for g in range(len(winners)):
            newstarts.append(starts[bedrecord][winners[g]])
            if len(winners) > 1:
                wins.append(str(starts[bedrecord][winners[g]]))
        newstart = int(round(sum(newstarts)/len(newstarts),0))
        newend = newstart + 1
        newscore = "_".join(wins)
        print >> newbed, chra + '\t' + str(newstart) + '\t' + str(newend) + '\t' + bedrecord + '\t' + newscore + '\t' + stranda
        k += 1
        continue
    coding = []
    count = -1
    for bios in biotypes[bedrecord]:
        count += 1
        if bios == 'protein_coding':
            coding.append(count)
    newstarts = []
    newscores = []
    for g in range(len(coding)):
    	newstarts.append(starts[bedrecord][coding[g]])
    	newscores.append(scores[bedrecord][coding[g]])
    winning = max(newscores)
    winners = []
    count = -1
    for score in newscores:
        count += 1
        if score == winning:
            winners.append(count)
    newerstarts = []
    wins = [str(winning)]
    for g in range(len(winners)):
        newerstarts.append(newstarts[winners[g]])
        if len(winners) > 1:
            wins.append(str(newstarts[winners[g]]))
    newstart = int(round(sum(newstarts)/len(newstarts),0))
    newend = newstart + 1
    newscore = "_".join(wins)
    print >> newbed, chra + '\t' + str(newstart) + '\t' + str(newend) + '\t' + bedrecord + '\t' + newscore + '\t' + stranda
    m += 1




    # if newbiotypes.count('protein_coding') == 0:
    # 	wins = []
    # 	wins.append(winning)
    # 	for newstart in newstarts:
    #     	wins.append(str(newstart))
    # 	newstarts = str(int(round(sum(newstarts)/len(newstarts),0)))
    # 	newends = str(int(round(sum(newends)/len(newends),0)))
    # 	newscores = '_'.join(wins)
    # 	print >> newbed, chra + '\t' + newstarts + '\t' + newends + '\t' + bedrecord + '\t' + newscores + '\t' + stranda
    #     m += 1
    # 	continue
    # if newbiotypes.count('protein_coding') == 1:
    # 	print >> newbed, chra + '\t' + str(newstarts[newbiotypes.index('protein_coding')]) + '\t' + str(newends[newbiotypes.index('protein_coding')]) + '\t' + bedrecord + '\t' + str(newscores[newbiotypes.index('protein_coding')]) + '\t' + strands[bedrecord][winner]
    #     p += 1
    #     continue

    # newerstarts = []
    # newerends = []
    # for d in range(len(coding)):
    # 	newerstarts.append(newstarts[coding])
    # 	newerends.append(newends[coding])
    # wins = []
    # wins.append(winning)
    # for start in newerstarts:
    #     wins.append(str(start))
    # newerstarts = str(int(round(sum(newerstarts)/len(newerstarts),0)))
    # newerends = str(int(round(sum(newerends)/len(newerends),0)))
    # newscores = '_'.join(wins)
    # print >> newbed, chra + '\t' + newerstarts + '\t' + newends + '\t' + bedrecord + '\t' + newscores + '\t' + stranda
    # r += 1

newbed.close()
print 'There were ' + str(j) + ' genes with 1 TSS.'
print 'There were ' + str(k) + ' genes with no protein coding TSS.'
print 'There were ' + str(m) + ' genes with at least 1 protein coding TSS.'

print 'Sorting records...'
sorter = 'sort -k1,1 -k2,2n -k3,3n -k4,4 ' + newbedfile + ' > ' + sortbed
sorting = subprocess.Popen(sorter,shell=True)
sorting.wait()