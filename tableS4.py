#!/usr/bin/env python
import glob

godir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/'
goreport = open('/mnt/lustre/home/cusanovich/Kd_Arrays/Supp/TableS4.txt','w')
reports = glob.glob(godir + '*_GOResults.txt')
genes = sorted([x.split('/')[-1].split('_')[0] for x in reports])
tablewidth = len(reports)

repodic = {}
godic = {}
for report in reports:
	currgene = report.split('/')[-1].split('_')[0]
	#print currgene
	currgo = open(report,'r')
	for line in currgo:
		if 'GO.ID' in line:
			continue
		liner = line.strip().split()

		if float(liner[-1]) >= 0.05:
			break
		gocater = line.strip().split('"')[3]
		if gocater not in repodic.keys():
			repodic[gocater] = ['NA']*tablewidth
			godic[gocater] = line.strip().split('"')[5]
		repodic[gocater][genes.index(currgene)] = "%.2f"%(float(liner[-4])/float(liner[-3])) + '(' + "%.2e"%float(liner[-1]) + ')'
	currgo.close()

print >> goreport, 'GO.ID\tGOTerm\t' + '\t'.join(genes)
for gocat in sorted(godic.keys()):
	print >> goreport, gocat + '\t"' + godic[gocat] + '"\t' + '\t'.join(repodic[gocat])

goreport.close()