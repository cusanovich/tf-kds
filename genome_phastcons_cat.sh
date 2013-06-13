#!/bin/bash
dir=/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/PhastCons/
outdir=`echo $dir | sed 's/PhastCons/FinalAnnots/'`
files=`ls ${dir}*.counts`

for file in $files
    do
        cat $file >> ${outdir}temp.phastcons.bed
    done

~/Programs/bin/sort-bed ${outdir}temp.phastcons.bed > ${outdir}master.phastcons.bed
rm ${outdir}temp.phastcons.bed