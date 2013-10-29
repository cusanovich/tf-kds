#!/usr/bin/env python
import os
import sys
import subprocess
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
from pybedtools import BedTool, featurefuncs
from kdfunc import pwmdicter

windowsize = 1000
windowname = str(windowsize/1000) + 'kb'
outdir = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/indiv_factor_beds/'
outboundonly = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_BoundOnly.bed'
oldout = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_old_sorted.bed'
newout = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_new_sorted.bed'
chiped = '/mnt/lustre/home/cusanovich/Kd_Arrays/EncodeChipSeq/sorted_EncodeChIP_Uniform_Combined.bed'
jacked = '/mnt/lustre/home/cusanovich/centipede/hg19_jack_centipede_sorted_pwms_clean.bed'
tsss = '/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed'
outfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_new.bed'
outsortfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/results' + windowname + '_combined_midpoint_new_sorted.bed'
newcombofile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/combinedChipandCenti_midpoint_new.bed'
newcombosortedfile = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/combinedChipandCenti_midpoint_new_sorted.bed'
pwms = open('/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt','r')
pwmlist = pwmdicter(pwms)
pwms.close()

def bedmaker(bedrecord,tssrecord,filerecord,unmatched=True):
	"""Writes bed-like record of all binding within 'w' kb of TSSs."""
	temper = open('./tempfile','w')
	sortrecord = bedrecord.sort()
	print >> temper, tssrecord.window(sortrecord, w=windowsize)
	if unmatched:
		print >> temper, tssrecord.window(sortrecord, w=windowsize, v=True)
		temper.close()
		sorter = 'grep . ./tempfile | sort -k1,1 -k2,2n -k3,3n -k4,4 > ' + filerecord + '; rm ./tempfile'
		sorting = subprocess.Popen(sorter,shell=True)
		sorting.wait()
	else:
		temper.close()
		sorter = 'grep . ./tempfile > ' + filerecord + '; rm ./tempfile'
		sorting = subprocess.Popen(sorter,shell=True)
		sorting.wait()
	return(0)

def newber(binner,factor):
	"""Returns bed record with all binding for a factor merged and renamed."""
	newbie = BedTool(binner).sort().merge(nms=True).each(featurefuncs.midpoint)
	newbie = newbie.each(featurefuncs.rename,factor)
	return(newbie)

print "Loading bed files..."
tssbed = BedTool(tsss)
jacker = BedTool(jacked)
chiper = BedTool(chiped)

print "Combining binding files..."
comber = chiper.cat(jacker,force_truncate=False,postmerge=False).moveto('./combtemp.bed')

#print "Generating oldstyle bed file..."
#oldstyle = comber.each(featurefuncs.midpoint)
#oldstyle2 = bedmaker(oldstyle,tssbed,oldout)

print "Generating newstyle bed files..."
for factor in sorted(pwmlist.keys()):
	print factor
	for model in pwmlist[factor]:
		print model
		grepper = 'grep ' + model + ' ./combtemp.bed >> ./temp1.bed'
		grepping = subprocess.Popen(grepper,shell=True)
		grepping.wait()
	newbie = newber('./temp1.bed',factor)
	newbie2 = bedmaker(newbie,tssbed,outdir + factor + '_midpoint_sorted.bed',unmatched=False)
	noob = newber('./temp1.bed',factor).moveto('./temp2.bed')
	cleaner = 'cat ./temp2.bed >> ' + newcombofile + '; rm ./temp1.bed; rm ./temp2.bed'
	cleaning = subprocess.Popen(cleaner,shell=True)
	cleaning.wait()

print "Combining newstyle bed files..."
beds = os.listdir(outdir)
for bed in beds:
	mover = 'cat ' + outdir + bed + ' >> ' + outfile
	moving = subprocess.Popen(mover,shell=True)
	moving.wait()

print "Cleaning up a bit..."
cleaner = 'sort -k1,1 -k2,2n -k3,3n -k4,4 ' + outfile + ' > ' + outsortfile + '; rm ' + outfile + '; rm ./combtemp.bed; sort -k1,1 -k2,2n -k3,3n -k4,4 ' + newcombofile + ' > ' + newcombosortedfile + '; rm ' + newcombofile
cleaning = subprocess.Popen(cleaner,shell=True)
cleaning.wait()