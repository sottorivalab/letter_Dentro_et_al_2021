source /home/salvatore.milite/.bashrc

conda activate phylowgs

#no copy number events in our simple simulation
rm -rf wgs_output

touch wgs_input/cnv_data.txt

mkdir wgs_output

python ../phylowgs/multievolve.py --num-chains 4 --ssms "./wgs_input/ssm_data.txt" --cnvs "./wgs_input/cnv_data.txt" --output-dir "./wgs_output" 