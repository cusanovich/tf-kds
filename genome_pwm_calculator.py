#!/usr/bin/env python
import sys
sys.path.insert(0,'/mnt/lustre/home/cusanovich/Programs/lib/python2.6/site-packages/pybedtools-0.6.2-py2.6-linux-x86_64.egg/pybedtools')
import pybedtools

d = pybedtools.BedTool("""
	chr1 840146 840165
	chr1 840733 840746""", from_string=True)


fasta = pybedtools.example_filename('/mnt/lustre/home/cusanovich/Genomes/allhg19_norandom.fasta')
d = d.sequence(fi=fasta)
d.seq("chr1:1-10",fasta)
d.print_sequence()