targets=`cat ../Annotations/TargetSummary.txt | awk -v FS="\t" -v OFS="," '{print $1,$2}'`
mkdir ../Results/RUV2_NSAveraged_Results
mkdir ../Outputs/RUV2_NSAveraged_Outputs
for target in $targets
    do
    factor=`echo $target | awk -v FS="," '{print $1}'`
    file=`echo $target | awk -v FS="," '{print $2}'`
    echo "Rscript KdLRTRUV2_AverageNS.R ${factor} ${file}" | qsub -l h_vmem=4g -o ../Outputs/RUV2_NSAveraged_Outputs/ -e ../Outputs/RUV2_NSAveraged_Outputs/ -wd `pwd` -N "KdAnalysis.${factor}"
done