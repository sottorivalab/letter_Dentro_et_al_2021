PIP_DIR="/Users/salvatore.milite/work/Trieste/analysis/letter_Dentro_et_al_2021/example_simulated/dpclust/inst/example/"

R --vanilla --slave -q -f "${PIP_DIR}dpclust_pipeline.R" --args -r 1 -d ./ dpclust_input/ -o ./ -i simulation_input.txt