#!/usr/bin/env python
import sys

if(len(sys.argv)!= 2):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: midpoints.py [Window Size]\n\n')
	sys.exit(1)

windowsize = sys.argv[1]
centifile = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'results_phastcons_sorted.bed','r')
outcentifile = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/' + windowsize + 'results_phastcons_sorted_midpoint.bed','w')

genes = []
for line in centifile:
	liner = line.strip().split()
	if len(liner) < 7:
		print >> outcentifile, "\t".join(liner)
		genes.append(liner[3])
		continue
	liner[7] = str(int(round((float(liner[7]) + float(liner[8]))/2,0)))
	liner[8] = str(int(liner[7]) + 1)
	middistance = int(liner[8]) - int(liner[2])
	if abs(middistance) < 10001:
		print >> outcentifile, "\t".join(liner)
		if liner[3] not in genes:
			genes.append(liner[3])
		continue
	if liner[3] not in genes:
		print >> outcentifile, "\t".join(liner[0:6])
		genes.append(liner[3])

centifile.close()
outcentifile.close()