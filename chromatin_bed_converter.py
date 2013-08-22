#!/usr/bin/env python
import sys
import subprocess

windowsize = '10000'
windowname = str(int(windowsize)/1000) + 'kb'

bindbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/combinedChipandCenti_midpoint_new_sorted.bed'
chrombed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/RawData/wgEncodeBroadHmmGm12878HMM.bed.gz'
openbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/RawData/wgEncodeAwgDnaseUwdukeGm12878UniPk.narrowPeak.gz'
opentagbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/DNase_tagged.bed'
chromcombobed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/DNase_Chromstates_combined_sorted.bed'
bindchrombed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/combinedbinding_chromstates.bed'
sortbindchrombed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/combinedbinding_chromstates_sorted.bed'
tssbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed'
allbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/allbinding_chromstatesandtss_overlap.bed'
sortallbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/' + windowname + '_allbinding_chromstatesandtss_overlap_sorted.bed'

tagger = "zcat " + openbed + " | sed 's/.\t0\t./0_Open/' > " + opentagbed + "; gzip " + opentagbed
cater = "zcat " + chrombed + " " + opentagbed + " | cut -f1-4 | sort -k1,1 -k2,2n -k3,3n -k4,4 > " + chromcombobed
binder = "cut -f1-5 " + bindbed + " | intersectBed -wa -wb -a stdin -b " + chromcombobed + " | cut -f1-9 > " + bindchrombed
sorter = "sort -k1,1 -k2,2n -k3,3n -k4,4 " + bindchrombed + " > " + sortbindchrombed
tsser = "windowBed -w " + windowsize + " -a " + tssbed + " -b " + sortbindchrombed + " > " + allbed + "; windowBed -v -w " + windowsize + " -a " + tssbed + " -b " + sortbindchrombed + " >> " + allbed
resorter = "sort -k1,1 -k2,2n -k3,3n -k4,4 " + allbed + " > " + sortallbed
cleanser = "rm " + bindchrombed + "; rm " + sortbindchrombed + "; rm " + allbed + "; rm " + opentagbed + ".gz; rm " + chromcombobed

print "Tagging DNase data..."
tagify = subprocess.Popen(tagger,shell=True)
tagify.wait()
print "Combine DNase and chromatin states..."
catify = subprocess.Popen(cater,shell=True)
catify.wait()
print "Overlapping data files..."
bindify = subprocess.Popen(binder,shell=True)
bindify.wait()
print "Sorting binding and chromatin states bed..."
sortify = subprocess.Popen(sorter,shell=True)
sortify.wait()
print "Finding binding near TSSs..."
tssify = subprocess.Popen(tsser,shell=True)
tssify.wait()
print "Re-sorting results..."
resortify = subprocess.Popen(resorter,shell=True)
resortify.wait()
print "Cleaning up after myself..."
cleansify = subprocess.Popen(cleanser,shell=True)
cleansify.wait()
print "Done."