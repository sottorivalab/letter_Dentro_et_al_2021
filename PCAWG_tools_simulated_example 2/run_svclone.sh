source /home/salvatore.milite/.bashrc

module load R

conda activate svclone

svclone cluster -s example --snv svclone_input/example_filtered_snvs.tsv --purity_ploidy svclone_input/pp_file.txt --out svclone_output -cfg ../svclone_config.ini