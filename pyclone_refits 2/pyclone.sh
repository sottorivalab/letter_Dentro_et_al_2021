#!/bin/zsh

# We assume to have a conda environment named pyclone-vi with all the dependencies already
# installed

# you should source your conda init script here
source ~/miniforge3/etc/profile.d/conda.sh

conda activate pyclone-vi

pyclone-vi fit -i pyclone_input.tsv -o pyclone_output.h5 -c 10 -d beta-binomial -r 25
pyclone-vi write-results-file -i pyclone_output.h5 -o pyclone_output.tsv
