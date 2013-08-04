#!/usr/bin/env python
import sys
import subprocess

windowsize = '5000'
windowname = '5kb'

bindbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/combinedChipandCenti_midpoint_sorted.bed'
chrombed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/RawData/wgEncodeBroadHmmGm12878HMM.bed.gz'
bindchrombed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/combinedbinding_chromstates.bed'
sortbindchrombed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/combinedbinding_chromstates_sorted.bed'
tssbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed'
allbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/allbinding_chromstatesandtss_overlap.bed'
sortallbed = '/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Overlaps/' + windowname + '_allbinding_chromstatesandtss_overlap_sorted.bed'

binder = "cut -f1-5 " + bindbed + " | intersectBed -wa -wb -a stdin -b " + chrombed + " | cut -f1-9 > " + bindchrombed
sorter = "sort -k1,1 -k2,2n -k3,3n -k4,4 " + bindchrombed + " > " + sortbindchrombed
tsser = "windowBed -w " + windowsize + " -a " + tssbed + " -b " + sortbindchrombed + " > " + allbed + "; windowBed -v -w " + windowsize + " -a " + tssbed + " -b " + sortbindchrombed + " >> " + allbed
resorter = "sort -k1,1 -k2,2n -k3,3n -k4,4 " + allbed + " > " + sortallbed
cleanser = "rm " + bindchrombed + "; rm " + sortbindchrombed + "; rm " + allbed

print "Overlapping data files..."
bindify = subprocess.Popen(binder,shell=True)
bindify.wait()
print "--> Binding and chromatin states overlapped..."
sortify = subprocess.Popen(sorter,shell=True)
sortify.wait()
print "--> Binding and chromatin states sorted..."
tssify = subprocess.Popen(tsser,shell=True)
tssify.wait()
print "---> TSSs found..."
resortify = subprocess.Popen(resorter,shell=True)
resortify.wait()
print "---> Overlaps sorted..."
cleansify = subprocess.Popen(cleanser,shell=True)
cleansify.wait()
print "Done."