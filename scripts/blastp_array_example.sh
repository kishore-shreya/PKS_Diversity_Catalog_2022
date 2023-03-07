#!/bin/bash
#SBATCH -J blastp
#SBATCH --output=blastp_%A_%a.out
#SBATCH --error=blastp_%A_%a.err
#SBATCH -p normal,hns
#SBATCH --time=500
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=20G

#!/bin/bash

cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs

module load biology
module load ncbi-blast+/

FILES=(*.fasta)
FILENAME=${FILES[${SLURM_ARRAY_TASK_ID}+999]}
#kept changing this number for set of 1000 tasks
# We have 18165 fasta files, so we have to run 18165 1-vs-all blastp queries. We do this using job arrays. # Each array can support up to 1000 jobs, so we have to submit 19 job arrays. # Task numbers in each array go from 0 to N, so we have to modify each slurm file, so that blastp is performed on different files in each array. The number we modify is the one in the line right above

blastp -db /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/blastdb/all_prot_alldbs_withClusterN_blast_db -query /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/${FILENAME} -outfmt 7 -num_threads 1 -out /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/blastp_scores/${FILENAME}.blastp.out -max_target_seqs 20000