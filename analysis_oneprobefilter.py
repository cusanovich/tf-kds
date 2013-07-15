#!/usr/bin/env python
import sys

if(len(sys.argv)!= 3):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: oneprobefilter.py [In File] [Out File]\n\n')
	sys.exit(1)

report = open(sys.argv[1],'r')
filtered = open(sys.argv[2],'w')
ensemblfile = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed','r')
ensemblgenes = set([x.strip().split()[3] for x in ensemblfile.readlines()])

x = 0
genes = {}
for line in report:
    if 'probeID' in line:
        continue
    gene = line.strip().split()[6]
    if gene not in ensemblgenes:
        x += 1
        continue
    if gene not in genes.keys():
        genes[gene] = line.strip().split()
        continue
    strand = line.strip().split()[5]
    if strand == '+' and genes[gene][1] > line.strip().split()[1]:
        continue
    if strand == '+' and genes[gene][1] < line.strip().split()[1]:
        genes[gene][1:4] = line.strip().split()[1:4]
    if strand == '-' and genes[gene][1] < line.strip().split()[1]:
        continue
    if strand == '-' and genes[gene][1] > line.strip().split()[1]:
        genes[gene][1:4] = line.strip().split()[1:4]
report.close()
print str(x) + ' genes were not in the CAGE dataset.'

print >> filtered,'chr\tstart\tend\tprobeID\tscore\tstrand\ttargetENSG\ttargetHGNC'
for gene in sorted(genes.iterkeys()):
    print >> filtered, '\t'.join(genes[gene])
filtered.close()
