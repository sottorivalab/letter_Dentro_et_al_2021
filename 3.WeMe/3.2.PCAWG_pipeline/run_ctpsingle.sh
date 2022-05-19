source /home/salvatore.milite/.bashrc

module load R

mkdir ctp_output

../CTPsingle/CTPsingle.R -f ./ctp_input/sample.snv.txt -o ./ctp_output/example_PCAWG_ctp -m ../CTPsingle/GammaAdjMatrices