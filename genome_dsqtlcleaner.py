#!/usr/bin/env python
import gzip

infile = gzip.open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/StartAnnots/GSE31388_eQtlTable.txt.gz','rb')
outfile = open('/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/StartAnnots/GSE31388_eQtlTable_cleaned.bed','w')
for line in infile:
	if '#' in line:
		continue
	liner = line.strip().split()
	print >> outfile, liner[0] + '\t' + liner[3] + '\t' + liner[3] + '\t' + liner[4] + ':' + liner[8] + '\t' + liner[5]
infile.close()
outfile.close()