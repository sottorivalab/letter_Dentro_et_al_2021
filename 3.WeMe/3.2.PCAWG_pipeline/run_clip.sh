source /home/salvatore.milite/.bashrc

conda activate clip

module load R

cd ../CliP

python3 run_clip_main.py -i $1 "../$1/clip_input/sample.snv.txt" "../$1/clip_input/sample.cna.txt" "../$1/clip_input/sample.purity.txt"