#!/bin/bash

dir=/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/
#outdir=`echo $dir | sed 's/PhastCons\///'`
files=`ls ${dir}*.bed`
for file in $files
    do
    	outfile=`echo $file | sed 's/.bed/.counts/'`
    	chra=`echo $file | sed 's/${dir}//'`
    	chrs=`echo $chra | sed 's/.phastCons46way.placental.bed//'`
        echo "~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/StartAnnots/hg19_jack_centipede_resorted.bed $file | grep -P ${chrs}\t > outfile" | qsub -l h_vmem=4g -o ../Outputs/ -e ../Outputs/ -wd `pwd` -N "phastcons_mean"
done


#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr1.phastCons46way.placental.bed | grep -P "chr1\t" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr1.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr2.phastCons46way.placental.bed | grep -P "chr2\t" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr2.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr3.phastCons46way.placental.bed | grep "chr3" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr3.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr4.phastCons46way.placental.bed | grep "chr4" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr4.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr5.phastCons46way.placental.bed | grep "chr5" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr5.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr6.phastCons46way.placental.bed | grep "chr6" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr6.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr7.phastCons46way.placental.bed | grep "chr7" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr7.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr8.phastCons46way.placental.bed | grep "chr8" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr8.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr9.phastCons46way.placental.bed | grep "chr9" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr9.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr10.phastCons46way.placental.bed | grep "chr10" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr10.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr11.phastCons46way.placental.bed | grep "chr11" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr11.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr12.phastCons46way.placental.bed | grep "chr12" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr12.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr13.phastCons46way.placental.bed | grep "chr13" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr13.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr14.phastCons46way.placental.bed | grep "chr14" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr14.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr15.phastCons46way.placental.bed | grep "chr15" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr15.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr16.phastCons46way.placental.bed | grep "chr16" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr16.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr17.phastCons46way.placental.bed | grep "chr17" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr17.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr18.phastCons46way.placental.bed | grep "chr18" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr18.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr19.phastCons46way.placental.bed | grep "chr19" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr19.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr20.phastCons46way.placental.bed | grep "chr20" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr20.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr21.phastCons46way.placental.bed | grep "chr21" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr21.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr22.phastCons46way.placental.bed | grep "chr22" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr22.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chrX.phastCons46way.placental.bed | grep "chrX" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chrX.phastCons46way.placental.counts;
#~/Programs/bin/bedmap --echo --mean --delim "\t" ~/Kd_Arrays/GenomeAnnotations/hg19_jack_centipede_resorted.bed ~/Kd_Arrays/GenomeAnnotations/PhastCons/chrY.phastCons46way.placental.bed | grep "chrY" > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chrY.phastCons46way.placental.counts;
