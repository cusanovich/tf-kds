#!/usr/bin/env python
import sys
import subprocess

if(len(sys.argv)!= 6):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: bedconverter.py [TSS BED File] [Centipede BED File] [ChIP BED File ] [Chromatin BED File] [Out Dir]\n\n')
	sys.exit(1)

tsss = sys.argv[1]
centi = sys.argv[2]
chip = sys.argv[3]
chrom = sys.argv[4]
olap = sys.argv[5]

#Pull appropriate fields from binding files and combine and sort for later steps
print "Cleaning binding files..."
centier = "zcat " + centi + " | cut -f1-5 - > " + olap + "all_centi.bed; cat " + olap + "all_centi.bed > " + olap + "allbinding_combined.bed"
chiper = "cut -f1-5 " + chip + " > " + olap + "all_chip.bed; cat " + olap + "all_chip.bed >> " + olap + "allbinding_combined.bed;"
sorter = "sort -k1,1 -k2,2n -k3,3n -k4,4 " + olap + "allbinding_combined.bed > " + olap + "allbinding_combined_sorted.bed"
centify = subprocess.Popen(centier,shell=True)
centify.wait()
print "--> Centipede done..."
chipify = subprocess.Popen(chiper,shell=True)
chipify.wait()
print "--> Encode done..."
sortify = subprocess.Popen(sorter,shell=True)
sortify.wait()
print "--> Binding sorted..."

print "Overlapping data files..."
binder = "zcat " + chrom + " | cut -f1-4 - | intersectBed -wa -wb -a " + olap + "allbinding_combined_sorted.bed -b stdin > " + olap + "allbinding_chromstates_overlap.bed"
tsser = "windowBed -w 10000 -a " + tsss + " -b " + olap + "allbinding_chromstates_overlap.bed > " + olap + "allbinding_chromstatesandtss_overlap.bed; windowBed -v -w 10000 -a " + tsss + " -b " + olap + "allbinding_chromstates_overlap.bed >> " + olap + "allbinding_chromstatesandtss_overlap.bed"
resorter = "sort -k1,1 -k2,2n -k3,3n -k4,4 " + olap + "allbinding_chromstatesandtss_overlap.bed > " + olap + "allbinding_chromstatesandtss_overlap_sorted.bed"
cleanser = "rm " + olap + "all_centi.bed; rm " + olap + "allbinding_combined.bed; rm " + olap + "all_chip.bed; rm " + olap + "allbinding_chromstates_overlap.bed; rm " + olap + "allbinding_chromstatesandtss_overlap.bed"
bindify = subprocess.Popen(binder,shell=True)
bindify.wait()
print "--> Binding and chromatin states overlapped..."
tssify = subprocess.Popen(tsser,shell=True)
tssify.wait()
print "---> TSSs found..."
resortify = subprocess.Popen(resorter,shell=True)
resortify.wait()
print "---> Overlaps sorted..."
cleansify = subprocess.Popen(cleanser,shell=True)
cleansify.wait()
print "Done."