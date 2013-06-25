#!/usr/bin/env python
import sys
import subprocess

indir = '/mnt/lustre/data/share/AllPwmSites/BcellAll/Results/Centipede/'
outfile = '/mnt/lustre/home/cusanovich/centipede/jack_centipede.bed'


motifs = []
lister = open('/mnt/lustre/home/cusanovich/centipede/src/AllUsablePWMs_AllBcell.txt','r')
for line in lister:
    motifs.append(line.strip())

lister.close()
motifs = {}.fromkeys(motifs).keys()
motifs.sort()

outter = open(outfile,'w')
for motif in motifs:
#    if motif[0:2] == 'MA':
#        continue

#    if motif[0:1] == 'S':
#        continue

    currbed = open(indir + motif + '.PostPr','r')
    for line in currbed:
        record = []
        line = line.strip().split()
        try:
            posterior = float(line[8])
        except ValueError:
            print "Bad line:"
            print motif
            print line
            continue
        if posterior < 0.95:
            continue

        record = line[0:6]
        #record.append(line[7])
        record[4] = line[8]
        print >> outter, "\t".join(record)

outter.close()

outdir = outfile.split('jack')[0]
sorter = 'sort -k1,1 -k2,2n -k3,3n -k4,4 ' + outfile + ' > ' + outdir + 'jack_centipede_sorted.bed; liftOver ' + outdir + 'jack_centipede_sorted.bed ~/liftover/hg18ToHg19.over.chain ' + outdir + 'hg19_jack_centipede_sorted.bed ' + outdir + 'jack_centipede_sorted_Unmapped.bed;'
sorting = subprocess.Popen(sorter,shell=True)
sorting.wait()

filtering = open(outdir + 'hg19_jack_centipede_sorted.bed','r')
filtered = open(outdir + 'hg19_jack_centipede_sorted_clean.bed','w')

previous = []
for line in filtering:
    liner = line.strip().split()
    if previous == []:
        previous = liner
        continue
    if liner[0:4] == previous[0:4]:
        previous[4] = str(max(float(previous[4]),float(liner[4])))
        previous[5] = "."
        continue
    print >> filtered, "\t".join(previous)
    previous = liner


#deleter =  'rm ' + outdir + 'jack_centipede.bed; rm ' + outdir + 'jack_centipede_sorted.bed; rm ' + outdir + 'hg19_jack_centipede_sorted.bed;'
#deleting = subprocess.Popen(deleter,shell=True)
#deleting.wait()