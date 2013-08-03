#!/usr/bin/env python
import sys

if(len(sys.argv)!= 3):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: resultsbed2factormatrix.py [BED File] [Results Matrix File]\n\n')
	sys.exit(1)

pwms = open('../Overlaps/allbinding_list.txt','r')
states = open('../Overlaps/chromstates_list.txt','r')
pather = '../Matrices/'
genes = open('../Overlaps/uniqueHT12genes.txt','r')

pwmlist = {}
for line in pwms:
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
genelist = [x.strip() for x in genelist]
genes.close()

for state in statelist:
    print state
    resultsbed = open(sys.argv[1],'r')
    resultname = pather + state + "." + str(sys.argv[2])
    resultsmatrix = open(resultname,'w')
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