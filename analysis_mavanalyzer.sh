targets=`cat /mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/TargetSummary.txt | egrep 'IRF4|HCST' |  awk -v FS="\t" -v OFS="," '{print $1,$2}'`

for target in $targets
    do
    factor=`echo $target | awk -v FS="," '{print $1}'`
    file=`echo $target | awk -v FS="," '{print $2}'`
    echo "Rscript /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/analysis_pub_mavolcanoes.R ${factor} ${file}" | qsub -l h_vmem=4g -o /mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Outputs/RUV2_NSAveraged_alt_Outputs/ -e /mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Outputs/RUV2_NSAveraged_alt_Outputs/ -wd `pwd` -N "PubMav.${factor}"
done