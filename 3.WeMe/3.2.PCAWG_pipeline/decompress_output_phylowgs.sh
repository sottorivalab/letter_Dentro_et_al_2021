source /home/salvatore.milite/.bashrc

conda activate phylowgs

cd wgs_output

python ../../phylowgs/write_results.py example ./trees.zip example.summ.json.gz example.muts.json.gz example_data.mutass.zip

cd ..