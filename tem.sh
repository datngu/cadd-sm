module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm


cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/cadd-sm

# env create -f environment.yml
conda activate cadd

bin=/mnt/ScratchProjects/Aqua-Faang/dat_projects/cadd-sm/bin

mkdir

wget ftp://ftp.ensembl.org/pub/release-66/emf/ensembl-compara/epo_6_primate/*

mkdir split_mod

./EMFSplitIndex.py -o split_mod/ *.gz

for i in {1..22} X
do 
    ./get_parameters_EPO.py -i whole_genome.fa -a split_mod -c $i -o $i.log
done
