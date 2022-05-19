source /home/salvatore.milite/.bashrc

conda activate clonehd

mkdir clonehd_output

cloneHD --snv  clonehd_input/clonehd_snvs.txt --pre clonehd_output/tumorSNV --seed 123 --trials 2\
 --nmax 5 --force --max-tcn 4 --purity clonehd_input/purity.txt --restarts 5 --mean-tcn clonehd_input/mean.tcn.txt