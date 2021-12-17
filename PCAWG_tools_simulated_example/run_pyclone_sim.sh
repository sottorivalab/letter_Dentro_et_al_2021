source ~/miniforge3/etc/profile.d/conda.sh

conda activate pyclone-vi


pyclone-vi fit -i simulation_input.tsv -o simulation.h5 -c 40 -d beta-binomial -r 10

pyclone-vi write-results-file -i simulation.h5 -o simulation.tsv