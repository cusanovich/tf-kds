targets=`ls /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Genotypes/ | awk -v FS="_" '{print $1}' | sort | uniq`
for target in $targets
    do
    echo "python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/eqtls_gemma_sidney.py ${target} DEandBound" | qsub -l h_vmem=2g -o /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Outputs/ -e /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Outputs/ -wd `pwd` -N "d${target}";
	echo "python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/eqtls_gemma_sidney.py ${target} Bound" | qsub -l h_vmem=2g -o /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Outputs/ -e /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Outputs/ -wd `pwd` -N "b${target}";
done