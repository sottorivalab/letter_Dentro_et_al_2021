source /home/salvatore.milite/miniconda3/etc/profile.d/conda.sh

conda activate pyclone-vi

mkdir output_pyclone

pyclone-vi fit -i simulation_input.tsv -o output_pyclone/simulation.h5 -c 40 -d beta-binomial -r 10

pyclone-vi write-results-file -i output_pyclone/simulation.h5 -o output_pyclone/simulation.tsv