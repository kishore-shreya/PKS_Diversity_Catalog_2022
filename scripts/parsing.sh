#!/bin/bash
#SBATCH --output=parsing_%j.out
#SBATCH --error=parsing_%j.err
#SBATCH -p normal,hns
#SBATCH --time=50
#SBATCH --nodes=1
#SBATCH -c 1

#!/bin/bash

cd /scratch/groups/khosla/Orphan_PKS/scripts/

perl /scratch/groups/khosla/Orphan_PKS/scripts/corrected_parse_blastp_sk.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/blastp_scores /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores 0 250 &

perl /scratch/groups/khosla/Orphan_PKS/scripts/corrected_parse_blastp_sk.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/blastp_scores /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores 250 500 &

perl /scratch/groups/khosla/Orphan_PKS/scripts/corrected_parse_blastp_sk.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/blastp_scores /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores 500 1000 &

#do it until index 18164

wait

exit 0
