#!/bin/bash
windowsize=1kb
window=1000
echo "Intersecting bed files..."
windowBed -a /mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed -b /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/master.phastcons.bed -w ${window} > /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/${windowsize}results_phastcons.bed;
windowBed -a /mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/HT12ensemblTSScombinedsorted.bed -b /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/master.phastcons.bed -w ${window} -v >> /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/${windowsize}results_phastcons.bed;
echo "Sorting bed files..."
sort -k1,1 -k2,2n -k3,3n -k4,4 /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/${windowsize}results_phastcons.bed > /mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/FinalAnnots/${windowsize}results_phastcons_sorted.bed;
echo "Calculating midpoints..."
python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/midpoints.py ${windowsize};