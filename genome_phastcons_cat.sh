#!/bin/bash
dir=/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/
outdir=`echo $dir | sed 's/PhastCons/FinalAnnots/'`
files=`ls ${dir}*.counts`

for file in $files
    do
        cat $file >> ${outdir}temp.phastcons.bed
    done

~/Programs/bin/sort-bed ${outdir}temp.phastcons.bed > ${outdir}master.phastcons.bed;
rm ${outdir}temp.phastcons.bed;
windowBed -a /mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12expr_ensemblTSScombinedsorted.bed -b /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/master.phastcons.bed -w 10000 > /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/10kbresults_phastcons.bed;
windowBed -a /mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12expr_ensemblTSScombinedsorted.bed -b /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/master.phastcons.bed -w 10000 -v >> /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/10kbresults_phastcons.bed;
sort -k1,1 -k2,2n -k3,3n -k4,4 /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/10kbresults_phastcons.bed > /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/10kbresults_phastcons_sorted.bed;
python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/midpoints.py;