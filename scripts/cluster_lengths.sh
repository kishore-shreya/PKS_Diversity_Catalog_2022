#!/bin/bash
#SBATCH --output=calculate_len.%j.out
#SBATCH --error=calculate_len.%j.err
#SBATCH -p normal,hns
#SBATCH --time=1000
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem-per-cpu=10G

#!/bin/bash

module load python/2.7.13
module load biology py-biopython/1.70_py27
module load math py-numpy
ml math py-scipy
ml math py-pandas
ml math py-sympy
ml viz py-matplotlib
module load python/2.7.13

python -u /scratch/groups/khosla/Orphan_PKS/scripts/cluster_length_sk.py -cluster_file /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.txt -inputf /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs -output_cluster_file /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt