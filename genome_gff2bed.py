#!/usr/bin/env python
import sys
import gzip
import subprocess

if(len(sys.argv)!= 4):
	sys.stderr.write('Error!! Supplied ' + str(len(sys.argv)) + ' Arguments!!\nUsage: gff2bed.py [Outdir] [In File] [Out File]\n\n')
	sys.exit(1)

gffin = gzip.open(sys.argv[1] + sys.argv[2],'rb')
bedder = sys.argv[1] + sys.argv[3]
bedout = open(bedder,'w')
#if '.bed.gz' not in sys.argv[2]:
#    bedout = gzip.open(sys.argv[2] + '.gz','w')
#else:
#    bedout = gzip.open(sys.argv[2],'w')
chosen = open('../qtlsources.txt','r')

choose = chosen.readlines()
choose = [x.strip() for x in choose]
chosen.close()

for line in gffin:
    liner = line.strip().split()
    if liner[1] not in choose:
        continue
    chrom = liner[0]
    start = str(int(liner[3]) - 1)
    end = liner[4]
    name = liner[2]
    score = liner[5]
    print >> bedout, chrom + '\t' + start + '\t' + end + '\t' + name + '\t' + score

gffin.close()
bedout.close()

outdir = sys.argv[1]
sorter = 'sort -k1,1 -k2,2n -k3,3n -k4,4 ' + bedder + ' > ' + outdir + 'sorted_' + sys.argv[3] + '; gzip ' + outdir + 'sorted_' + sys.argv[3]
subprocess.Popen(sorter,shell=True)
