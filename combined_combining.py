#!/usr/bin/env python
import os
import sys
import subprocess
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
import pybedtools

windowsize = 10000
windowname = str(windowsize/1000) + 'kb'
outter = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/combinedChipandCenti_midpoint.bed'
sortcombo = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/combinedChipandCenti_midpoint_sorted.bed'
outbounder = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint.bed'
outboundonly = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_BoundOnly.bed'
finalsortcombo = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_sorted.bed'
chiper = '/mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/sorted_EncodeChIP_Uniform_Combined.bed'
jacked = '/mnt/lustre/home/cusanovich/centipede/hg19_jack_centipede_sorted_pwms_clean.bed'
tsss = '/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed'


if 'combinedChipandCenti_midpoint_sorted.bed' not in os.listdir('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/'):
	print 'Calculating midpoints and combining files...'
	outcombo = open(outter,'w')
	jacks = open(jacked,'r')
	for line in jacks:
		if line == "":
			continue
		liner = line.strip().split()
		liner[1] = str(int(round((float(liner[1]) + float(liner[2]))/2,0)))
		liner[2] = str(int(liner[1]) + 1)
		print >> outcombo, "\t".join(liner)

	jacks.close()

	chips = open(chiper,'r')
	for line in chips:
		if line == "":
			continue
		liner = line.strip().split()
		liner[1] = str(int(round((float(liner[1]) + float(liner[2]))/2,0)))
		liner[2] = str(int(liner[1]) + 1)
		print >> outcombo, "\t".join(liner)

	chips.close()
	outcombo.close()

	print 'Sorting combined records...'
	sorter = 'grep . ' + outter + ' | sort -k1,1 -k2,2n -k3,3n -k4,4 > ' + sortcombo + '; rm ' + outter
	sorting = subprocess.Popen(sorter,shell=True)
	sorting.wait()

print 'Intersecting records with TSSs...'
outboundbed = open(outbounder,'w')
outjustboundbed = open(outboundonly,'w')
tssbed = pybedtools.BedTool(tsss)
combobed = pybedtools.BedTool(sortcombo)

print >> outboundbed, tssbed.window(combobed, w=windowsize)
print >> outjustboundbed, tssbed.window(combobed, w=windowsize)
print >> outboundbed, tssbed.window(combobed, w=windowsize, v=True)

outboundbed.close()
outjustboundbed.close()

print 'Sorting intersected records...'
sorter = 'grep . ' + outbounder + ' | sort -k1,1 -k2,2n -k3,3n -k4,4 > ' + finalsortcombo + '; rm ' + outbounder
sorting = subprocess.Popen(sorter,shell=True)
sorting.wait()