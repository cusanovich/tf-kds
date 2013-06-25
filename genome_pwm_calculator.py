#!/usr/bin/env python
import sys
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
import pybedtools
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/Bio')
import Bio
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/Bio')
from Bio import motifs
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy
from datetime import datetime
fasta = pybedtools.BedTool('/mnt/lustre/data/share/HumanGenome/allhg18_norandom.fasta')
matrices = open('/data/share/TRANSFAC/2011.3_nb/dat/matrix.dat','r')
bed = open('/mnt/lustre/home/cusanovich/centipede/jack_centipede_sorted.bed','r')
outbed = open('/mnt/lustre/home/cusanovich/centipede/jack_centipede_sorted_pwms_timedright2.bed','w')
#liner = ['chr1','847521','847534','M01066','0.9675','-']

d = pybedtools.BedTool("""chr1 840146 840165""", from_string=True)
genome = d.sequence(fi=fasta)
test = motifs.parse(matrices,"TRANSFAC")
motifers = []
for i in range(len(test)):
	test[i].pseudocounts = 0.5
	motifers.append(test[i]['AC'])

matrices.close()
x=1
pssms = {}
for line in bed:
	liner = line.strip().split()
	motifed = motifers.index(liner[3])
	#For some reason, Jack and Roger's bed files seem to be off in coordinates?!?!?!
	testbed = liner[0] + ':' + str(int(liner[1])-2) + '-' + str(int(liner[2])+2)
	testseq = Seq(genome.seq(testbed,fasta),IUPAC.unambiguous_dna)
	if liner[5] == '-':
		testseq = testseq.reverse_complement()

	if liner[3] in pssms.keys():
		pwms = pssms[liner[3]].calculate(testseq)
	else:
		pwms = test[motifed].pssm.calculate(testseq)
		pssms[liner[3]] = test[motifed].pssm
	
	maxscore = max(pwms)
	maxind = numpy.argmax(pwms) - 2
	pwmlen = len(test[motifed].pwm[1]) - 1
	if liner[5] != '-':
		newstart = str(int(liner[1]) + maxind)
		newend = str(int(newstart) + pwmlen)
	if liner[5] == '-':
		newend = str(int(liner[2]) - maxind)
		newstart = str(int(newend) - pwmlen)

	liner[1] = newstart
	liner[2] = newend
	print >> outbed, "\t".join(liner) + "\t" + str(maxscore)
	if x%10000 == 0:
		print str(datetime.now()) + " --> %s records processed..." % x

	x+=1

bed.close()
outbed.close()
