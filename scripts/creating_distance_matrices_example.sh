#!/bin/bash
#SBATCH --output=matrix6_%j.out
#SBATCH --error=matrix6_%j.err
#SBATCH -p normal,hns
#SBATCH --time=2880
#SBATCH --nodes=1
#SBATCH -c 2

#!/bin/bash

perl /scratch/groups/khosla/Orphan_PKS/scripts/creating_distance_matrices_by_year.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ClusterTablesTimeline/clusterTable2019.nonidentical.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/DistanceMatrices/distanceMatrix2019.txt &

perl /scratch/groups/khosla/Orphan_PKS/scripts/creating_distance_matrices_by_year.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ClusterTablesTimeline/clusterTable2020.nonidentical.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/DistanceMatrices/distanceMatrix2020.txt &

wait

exit 0

#did this from 1994 all the way to 2022