#!/usr/bin/env python
import os
import sys
import subprocess
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
from pybedtools import BedTool, featurefuncs

windowsize = 10000
windowname = str(windowsize/1000) + 'kb'
indir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/StartAnnots/'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/'
annots = ['GSE31388_eQtlTable_cleaned.bed','sorted_PritchardQTLs_merged.bed','gwascatalog_ucsc_merged.bed']
#annots = ['GSE31388_eQtlTable_cleaned.bed']
jacked = BedTool('/mnt/lustre/home/cusanovich/centipede/hg19_jack_centipede_sorted_pwms_clean.bed')
tss = BedTool('/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed')

for annot in annots:
	print annot
	currannot = BedTool(indir + annot)
	currout = annot.split('.')[0]
	print 'Intersecting...'
	inter = jacked.intersect(currannot,wa=True,wb=True).moveto(outdir + windowname + '_' + currout + '_centipede_intersect.bed')
	inter = BedTool(outdir + windowname + '_' + currout + '_centipede_intersect.bed')
	print 'Calculating midpoints...'
	intermid = inter.each(featurefuncs.midpoint).moveto(outdir + windowname + '_' + currout + '_centipede_intersect_midpoint.bed')
	print 'Finding TSSs...'
	inter = BedTool(outdir + windowname + '_' + currout + '_centipede_intersect_midpoint.bed')
	outter = tss.window(intermid,w=windowsize).moveto(outdir + windowname + '_' + currout + '_insite.bed')