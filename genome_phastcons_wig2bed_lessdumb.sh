#!/bin/bash

dir=/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/
#outdir=`echo $dir | sed 's/PhastCons\///'`
files=`ls ${dir}*.wigFix.gz`
for file in $files
    do
    	outfile=`echo $file | sed 's/.wigFix.gz/.bed/'`
        echo "zcat $file | /mnt/lustre/home/cusanovich/Programs/bin/wig2bed - | /mnt/lustre/home/cusanovich/Programs/bin/sort-bed - > $outfile" | qsub -l h_vmem=4g -o ../Outputs/ -e ../Outputs/ -wd `pwd` -N "wig2bed"
done



# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr1.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr1.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr2.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr2.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr3.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr3.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr4.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr4.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr5.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr5.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr6.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr6.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr7.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr7.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr8.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr8.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr9.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr9.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr10.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr10.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr11.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr11.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr12.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr12.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr13.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr13.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr14.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr14.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr15.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr15.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr16.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr16.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr17.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr17.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr18.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr18.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr19.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr19.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr20.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr20.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chr22.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chr22.phastCons46way.placental.bed;
# zcat /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/chrX.phastCons46way.placental.wigFix.gz | wig2bed - | sort-bed - > ~/Kd_Arrays/GenomeAnnotations/PhastCons/chrX.phastCons46way.placental.bed;