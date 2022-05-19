source /home/salvatore.milite/.bashrc

mkdir sclust_output


cd sclust_input

../../Sclust/bin/Sclust cluster -i example

cd ..

mv sclust_input/example_mclusters.txt sclust_input/example_cluster_assignments.txt sclust_output