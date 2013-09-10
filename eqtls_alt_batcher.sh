targets=`ls /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Genotypes/ | awk -v FS="_" '{print $1}' | sort | uniq`
for target in $targets
    do
    echo "python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/eqtls_alt_gemma_sidney.py ${target} DEandBound" | qsub -l h_vmem=2g -o /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltOutputs/ -e /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltOutputs/ -wd `pwd` -N "ad${target}";
	echo "python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/eqtls_alt_gemma_sidney.py ${target} Bound" | qsub -l h_vmem=2g -o /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltOutputs/ -e /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltOutputs/ -wd `pwd` -N "ab${target}";
	echo "python /mnt/lustre/home/cusanovich/Kd_Arrays/Scripts/eqtls_alt_gemma_sidney.py ${target} Unbound" | qsub -l h_vmem=2g -o /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltOutputs/ -e /mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/AltOutputs/ -wd `pwd` -N "au${target}";
done