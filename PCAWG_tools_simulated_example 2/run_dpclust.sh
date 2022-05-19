source /home/salvatore.milite/.bashrc



PIP_DIR="/group/sottoriva/salvatore_milite/letter_Dentro_et_al_2021/PCAWG_tools_simulated_example/dpclust/inst/example/"

module load R


R --vanilla --slave -q -f "${PIP_DIR}dpclust_pipeline.R" --args -r 1 -d ./ -o ./ -i simulation.txt