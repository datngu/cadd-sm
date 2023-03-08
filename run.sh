#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=Brain   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm






run_nextflow () {
  tis=$1
  genome=/mnt/users/ngda/genomes/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa
  annotation=/mnt/users/ngda/genomes/atlantic_salmon/Salmo_salar.Ssal_v3.1.106.gtf
  #meta
  meta_file=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/meta/${tis}.csv
  # genotype
  genotype=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/WGS/genotype.vcf.gz
  # atac
  atac_bam=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/${tis}/atac_bam/*{.bam,.bai}
  atac_count=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/${tis}/atac_consensus_peak_featureCounts.txt
  # rna
  rna_bam=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/${tis}/rna_bam/*{.bam,.bai}
  rna_count=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/${tis}/rna_gene_level_count_salmon.txt    
  # nextflow configs
  nextflow_res_dir=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/results/${tis}
  nextflow_work_dir=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/nextflow_work/${tis}
  mkdir $nextflow_work_dir
  cd $nextflow_work_dir
  git clone https://github.com/datngu/nf-rasqual-dev.git
  cd nf-rasqual-dev
  git pull
  export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity
  nextflow run main.nf -resume --meta $meta_file --genome $genome --annotation $annotation --atac_bam $atac_bam --atac_count $atac_count --rna_bam $rna_bam --rna_count $rna_count --genotype $genotype --outdir $nextflow_res_dir
}

#tissue_list="Brain Gill Gonad Liver Muscle"
run_nextflow Brain