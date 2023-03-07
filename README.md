# PKS_Diversity_Catalog_2022

########################################################################################
#########################################################################################
 
                                 DISTINCT NON-REDUNDANT PKS CATALOG (2022)

########################################################################################
#########################################################################################

**Note all mentioned scripts below are in the directory labelled scripts above**


The updated disitnct non-redundant polyketide synthases (PKS) catalog and further analysis was performed for the review article
"Diversity of Assembly Line Polyketide Synthases" written by Shreya Kishore (a,b) and Chaitan Khosla (a,b,c) <br />
a Department of Chemistry; b Stanford ChEM-H; c Department of Chemical Engineering; Stanford University, Stanford, CA 94305, USA in 2023.

More details can be found here: <>

This work was largely based on the one performed for the original article "Computational identification and analysis of orphan assembly-line polyketide synthases" written by Robert V O’Brien (1), Ronald W Davis (2,3), Chaitan Khosla (1,2,4) and Maureen E Hillenmeyer (2,3)<br /> 
1. Department of Chemistry, Stanford University, Stanford, CA, USA; 
2. Department of Biochemistry, Stanford University, Stanford, CA, USA; 
3. Stanford Genome Technology Center, 855 South California Avenue, Palo Alto, CA, USA; and 
4. Chemical Engineering, Stanford University, Stanford, CA, USA
published in The Journal of Antibiotics in 2014. As well as on the updated analysis as shown in the follwoing article "Evolution and Diversity of Assembly-Line Polyketide Synthases" written by Aleksandra Nivina (a,b), Kai P. Yuet (a,b), Jake Hsu (b,c) and Chaitan Khosla (a,b,c) <br />
a Department of Chemistry; b Stanford ChEM-H; c Department of Chemical Engineering; Stanford University, Stanford, CA 94305, USA for Chemical Reviews journal, in 2019.

Copyright (C) 2022 Maureen Hillenmeyer, Jake Hsu, Aleksandra Nivina, Shreya Kishore

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

To see the GNU General Public License, Plesase see 
<http://www.gnu.org/licenses/>.

Part of the work was performed on Stanford Sherlock Cluster. The rest was performed locally.

#########################################################################################
#########################################################################################
The first part of the analysis was performed using scripts from 2013.

##########################################################################################################################
*****

Install miniconda
Downloaded NCBI blast databases August 4th, 2022
https://ftp.ncbi.nih.gov/blast/db/
2.7.1

Blast consensus sequences against these dbs
1. Env_nt (DNA, environmental samples)
2. Nt (DNA, Partially non-redundant nucleotide sequences from all traditional divisions of GenBank, EMBL, and DDBJ excluding GSS,STS, PAT, EST,                HTG, and WGS)
3. Patnt (DNA, Nucleotide sequences derived from the Patent division of GenBank)
4. Tsa_nt (DNA, Transcriptome Shotgun Assembly (TSA) sequences)
5. ref_prok_rep_genomes (DNA, Refseq prokaryote representative genomes (contains refseq assembly)) 
6. ref_euk_rep_genomes (DNA, RefSeq Eukaryotic Representative Genome Database)
7. refseq_genomic (v4 db: until Feb 2020) 
8. other_genomic (v4 db:until Feb 2020)

WGS is no longer a blast db, SRA section explained later

### Running tblastn to Find KSs for the 8 BLAST Databases
```
module load biology 
module load ncbi-blast 
tblastn -query /scratch/groups/khosla/Orphan_PKS/KSSignatureConsensusPKSDB.fasta -db /scratch/groups/khosla/Orphan_PKS/blastdb/tsa_nt/tsa_nt -out KSSignature_tsa_nt.blastout.100000alignments.evalue100 -num_alignments 100000 -seg no -evalue 100
```

Repeat for the remaining BLAST databases

### Workaround for WGS: 
For wgs need to use the workaround because there's just too much data on wgs to search all wgs since 2016 
Workaround: ftp://ftp.hgc.jp/pub/mirror/ncbi/blast/WGS_TOOLS/README_BLASTWGS.txt <br />
So basically, taxid2wgs gets an alias file of a list acessions that is recongized by tblastn_vdb
tblastn_vdb it takes the list of accessions and tries to fetch the sra files from NCBI on the fly before running blast

TaxID of Organisms that I am Interested in: <br />
- Archae: 2157 <br />
- Bacteria: 2 <br />
- Eukaryota: 2759 <br />

To get the list of accessions: <br />
```
ml perl 
cpanm install LWP::Protocol::https 
cpanm install LWP::UserAgent 
perl /scratch/groups/khosla/Orphan_PKS/blastdb/WGS/taxid2wgs.pl -title "Archaea WGS" -alias_file archaea-wgs 2157
```

Created a file called: archaea-wgs.nvl <br />
Repeated for bacteria and eukaryota <br />

### tblastn_vdb to Find KSs for the SRA WGS Sequences

```
/scratch/groups/khosla/Orphan_PKS/ncbi-blast-2.13.0+/bin/tblastn_vdb -query /scratch/groups/khosla/Orphan_PKS/KSSignatureConsensusPKSDB.fasta -db /scratch/groups/khosla/Orphan_PKS/blastdb/WGS/x00 -out x00_KSSignature_wgs_archaea.blastout.100000alignments.evalue100 -num_alignments 100000 -seg no -evalue 100
```

### Cluster KS
force them to have at least 3kb between HSPs to count as a unique one. 
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/clusterKSs.pl /scratch/groups/khosla/Orphan_PKS/blastdb/env_nt/KSSignature_env_nt.blastout.100000alignments.evalue100 >KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart
```
Repeat for the remaining BLAST databases

######################################################################################################
# *****
### Filter genomes with at least 3 KSs
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/filterGenomes.pl KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart 3 > KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3 
```
Repeat for the remaining BLAST databases and SRA WGS

### How many entries meet the criteria?

cut -f 1 KSSignature_clusters.patnt.max.20000.min.3000.bp.apart.min3 |sort -u |wc <br />
383  383 4191 <br />

Repeat for the remaining BLAST databases and SRA WGS

| | | | |
| :---| :----|  :--- |:----   | 
| 196 |  196 |  3091 |				1. Env_nt | 
| 4341 |  4341 |  47928 | 				 2. Nt |
|383    |383  |  4191	|				   3. Patnt|
0   |   0      |0							 |   4. Tsa_nt |
254  |  254   | 3459					|   5. ref_euk_rep_genomes |
3120  | 3120  | 55250				|	  6. ref_prok_rep_genomes |
27852  |27852|  479501		  |  7. refseq_genomic (v4 db: until Feb 2020) |
4347 |  4347  | 65614 			|	  8. other_genomic (v4 db:until Feb 2020)|
1438  | 1438   |20132				 |  9. refseq_genomic_sinceFeb2020 |
1	   |  1		 |   18					    |10. WGS_archaea|
3398  | 3398 |  56703				  |11.  WGS_eukaryota|
898765 |898765| 15707662			|12. WGS_bacteria|

**Total NO. Of NCBI records with >=3 KS domains clustered within 20kb = 978,883** <br />
**Total NO. Of PKS sequences in the 978,883 records = 938,028**

# ****************************************************************************************************** 
### Fetch Genbanks
```
ml perl
cpanm install Bio::DB::EUtilities
perl /scratch/groups/khosla/Orphan_PKS/scripts/fetchGenbanks.pl KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3 gbwithparts
```
Num of ids = 196 <br />
[file called KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3.gbwithparts is created]

Repeat for the remaining BLAST databases and SRA WGS

### Separate Genbank File into Individual Files
```
ml perl
cpanm install Bio::SeqIO
perl /scratch/groups/khosla/Orphan_PKS/scripts/separateSequenceFiles.pl /scratch/groups/khosla/Orphan_PKS/KS_Cluster_results/KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3.gbwithparts genbank genbank.env_nt.min3.newset1
```
creates new folder called "genbank.env_nt.min3" and puts all individual files in there <br />
Repeat for the remaining BLAST databases and SRA WGS

### Convert Genbanks into FASTA for antiSMASH

Need to run this in the genbank.env_nt.min3.newset1 directory
```
ml perl
cpanm install Bio::SeqIO
for i in *.gb; do
perl /scratch/groups/khosla/Orphan_PKS/scripts/gb2fasta.pl /scratch/groups/khosla/Orphan_PKS/KS_Cluster_results/genbank.env_nt.min3.newset1/$i /scratch/groups/khosla/Orphan_PKS/KS_Cluster_results/env_nt_fasta/${i%%.*}.fasta
done
```
Repeat for the remaining BLAST databases and SRA WGS

##########################################################################################################################
# *****

### Install AntiSMASH Locally on Sherlock

Here is description of how to locally download antiSMASH: https://docs.antismash.secondarymetabolites.org/install/ 

1. Download anaconda first: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html 
2. Tried installing fully via bioconda but did not work. Too many issues. So the only way it worked was to install dependencies via bioconda and then pip install antismash. When I installed dpendencies via bioconda, and tried to download the antismash databases, it gave me an error on diamond package. So I had to:
```conda install hmmer2 hmmer diamond==0.8.36 fasttree prodigal blast muscle==3.8.1551 glimmerhmm meme==4.11.2 ```
```conda install -c conda-forge openjdk ```

### Run antiSMASH
```
for i in *; do
	antismash $i --genefinding-tool prodigal --output-dir /scratch/khosla/groups/Orphan_PKS/results/antismash/env_nt/${i%%.*} 
done
 ```
############################################################################################################

### Build description file and Get dates for the genbank files
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/build_accession_desc_jh.pl /scratch/groups/khosla/Orphan_PKS/results/KS_Cluster_results_min3/KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3 > /scratch/groups/khosla/Orphan_PKS/results/desc_files/KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3.accession_desc.txt
 ```
 Repeat for the remaining BLAST databases and SRA WGS
 
 catenated all into KSSignature_clusters.alldbs.max.20000.min.3000.bp.apart.min3.accession_desc.txt
 ```
perl /scratch/groups/khosla/Orphan_PKS/scripts/getSequenceDates.pl /scratch/groups/khosla/Orphan_PKS/results/genbank_files/env_nt/genbank.env_nt.min3.newset1 > /scratch/groups/khosla/Orphan_PKS/results/accession_dates /accession_dates_env_nt.txt
 ```
 Repeat for the remaining BLAST databases and SRA WGS
 
 catenated all into accession_dates_alldbs.txt
 
 
### Extract interesting clusters from all databases; here env_nt shown as example.

Repeat for the remaining BLAST databases and SRA WGS

**1. Get interesting clusters (at least 3 KS)**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/extractInterestingAntismashClusters_sk.pl /scratch/groups/khosla/Orphan_PKS/results/desc_and_accessionDate_files/env_nt/KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3.accession_desc.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/env_nt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/env_nt 3 > /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/env_nt/extractInterestingAntismashClusters_env_nt_v6.output
```
**2. Get cluster sequences**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/getClusterSequence.pl /scratch/groups/khosla/Orphan_PKS/results/genbank_files/env_nt/genbank.env_nt.min3.newset1 /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interestingClusters.env_nt.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt > /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interesting_cluster_seq_env_nt.txt
```
Combine all interesting_cluster_seq_*.txt > interesting_cluster_seq_alldbs.txt

**3. Find synonyms**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/findRedundantClusters_sk_parallel.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interesting_cluster_seq_env_nt.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt > /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interestingClusterSynonyms.env_nt.txt
```

Note: modified script to findRedundantClusters_sk_parallel to work in half the time because otherwise the predicted time was ~6000 days which is nearly impoosible. New script will compare each set of two only once instead of twice like it used to before. Can also take input parameter to process subsets of the data but compare it against the entire set. This way I can run multiples in parallel and then cat all the files together in the end.

**4. Remove redundant**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/removeRedundantClusters.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interestingClusters.env_nt.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interestingClusterSynonyms.env_nt.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt /scratch/groups/khosla/Orphan_PKS/results/desc_files/env_nt/KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3.accession_desc.txt /scratch/groups/khosla/Orphan_PKS/results/accession_dates/accession_dates_env_nt.txt > /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/interestingClustersNonRedundant.env_nt.txt 
```

**How many non-redundant PKS Clusters so far from each of the databases?**
| | |
| :--- |:----   | 
|160		|1. Env_nt|
|6564	|2. Nt|
|174		|3. Patnt|
|0		|4. Tsa_nt|
|47		|5. ref_euk_rep_genomes|
|3188	|6. ref_prok_rep_genomes |
|24036	 |7. refseq_genomic (v4 db: until Feb 2020) |
|	5037	|8. other_genomic (v4 db:until Feb 2020)|
|	2393	|9. refseq_genomic_sinceFeb2020 |
|	1		|10. WGS_archaea|
|	588	        |11.  WGS_eukaryota|
|	3554       |12. WGS_bacteria|
|45742	  |All together!!|

Note: for some of them like ref_prok_rep_genomes and refseq_genomic_sinceFeb2020 the number of interesting non-redundant clusters found is more than the number after 3 KS clustering before antiSMASH. This is because antiSMASH on a single ncbi accession can pull more than one interesting PKS cluster.

**5. Make summary cluster table (nonredundant)**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/makeClusterTableSummary_sk.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/interestingClustersNonRedundant.alldbs.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/alldbs /scratch/groups/khosla/Orphan_PKS/results/genbank_files/alldbs /scratch/groups/khosla/Orphan_PKS/results/desc_and_accessionDate_files/accession_dates_alldbs.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/interesting_cluster_seq_alldbs.txt > /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundant.alldbs.txt
```

**6. Consolidate similar**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/identifySimilarClusters.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundant.alldbs.txt > /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/SimilarClusters.alldbs.txt
```

**7. Cluster table nonsimilar non redundant**
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/makeClusterTableSummaryNonSimilar.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundant.alldbs.txt /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/SimilarClusters.alldbs.txt /scratch/groups/khosla/Orphan_PKS/results/desc_and_accessionDate_files/KSSignature_clusters.alldbs.max.20000.min.3000.bp.apart.min3.accession_desc.txt /scratch/groups/khosla/Orphan_PKS/results/desc_and_accessionDate_files/accession_dates_alldbs.txt > /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.txtfindRedundantClusters.
```

#########################################################################################
#########################################################################################
# Starting from here, the workflow uses scripts from 2018.

**8. Extracting biosynthetic enzyme protein sequences from GenBank to fasta**
use python 2.7 <br />
downgrade biopython <br />
bash script:
```
module load python/2.7.13
module load biology py-biopython/1.70_py27
module load math py-numpy
ml math py-scipy
ml math py-pandas
ml math py-sympy
ml viz py-matplotlib
module load python/2.7.13

python -u /scratch/groups/khosla/Orphan_PKS/scripts/extracting_protein_seq_sk.py -cluster_file /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt/interesting_clusters/clusterTableNonRedundantNonSimilar.env_nt.txt -startID 0 -endID 60 -inputf /scratch/groups/khosla/Orphan_PKS/results/antismash_results_v4/env_nt -outputf /scratch/groups/khosla/Orphan_PKS/results/prot_fasta_withClusterN_for_blastp/env_nt_v4 &

wait

exit 0
```
**9. Generating a GenBank database to perform 1-vs-all blastp queries**
Generating a GenBank database to perform 1-vs-all blastp queries.

```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs 
cat *.fasta > allfiles.fasta
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs
mkdir blastdb

makeblastdb -in /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/allfiles.fasta -dbtype prot -out /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/blastdb/all_prot_alldbs_withClusterN_blast_db
mv allfiles.fasta output/ 
```
**10. Running blastp**
We have 18164 fasta files, so we have to run 18164 1-vs-all blastp queries. We do this using job arrays. Each array can support up to 1000 jobs, so we have to submit 19 job arrays. Task numbers in each array go from 0 to N, so we have to modify each slurm file, so that blastp is performed on different files in each array. A user can only submit 1 job array of 1000 tasks, so I submitted them sequentially.

Example sbatch script is located in the scripts folder called blastp_array_example.sh 

sequentially run:
sbatch --array=0-999 blastp_array.sh <br />
sbatch --array=0-999 blastp_array_999.sh <br />
sbatch --array=0-999 blastp_array_1998.sh <br />
sbatch --array=0-999 blastp_array_2997.sh <br />
sbatch --array=0-999 blastp_array_3996.sh <br />
sbatch --array=0-999 blastp_array_4995.sh <br />
sbatch --array=0-999 blastp_array_5994.sh <br />
sbatch --array=0-999 blastp_array_6993.sh <br />
sbatch --array=0-999 blastp_array_7992.sh <br />
sbatch --array=0-999 blastp_array_8991.sh <br />
sbatch --array=0-999 blastp_array_9990.sh <br />
sbatch --array=0-999 blastp_array_10989.sh <br />
sbatch --array=0-999 blastp_array_11988.sh <br />
sbatch --array=0-999 blastp_array_12987.sh <br />
sbatch --array=0-999 blastp_array_13986.sh <br />
sbatch --array=0-999 blastp_array_14985.sh <br />
sbatch --array=0-999 blastp_array_15984.sh <br />
sbatch --array=0-999 blastp_array_16983.sh <br />
sbatch --array=0-185 blastp_array_17982.sh <br />
(total 18164)

**11. Calculating query lengths**
Because we now have to calculate everything per amino acid length and not per nucleotide length, we have to modify the clusterTable to reflect the amino acid length.
```
sbatch cluster_lengths.sh
```

Resulting table goes to /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt

**12.  Parsing blastp output files**
```
sbatch parsing.sh
```

**13. Finding redundant clusters and chaining them**
This script finds redundant clusters based on their similarity scores

```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/
perl /scratch/groups/khosla/Orphan_PKS/scripts/findRedundant_and_chain_blastp.pl parsed_scores_alldbs 0.9
```

I got: <br />
non-redundant = 6329 <br />
redundant sets = 2470 <br />

### Total: 8799 non redundant, unique PKS clusters!!

**14. Getting the correct domain numbers**
It turns out that the original script mis-counted the number of C domains: it took these numbers from svg file, where C domain and C-terminal docking domains are labeled in a very similar way. Aleks wrote a script that gets the numbers of domains from GenBank files, where this information is correct as annotated by antiSMASH. I modified this script to work with new antiSMASH 6 annotation.

```
module load python/2.7.13
module load biology py-biopython/1.70_py27
module load math py-numpy
ml math py-scipy
ml math py-pandas
ml math py-sympy
ml viz py-matplotlib
module load python/2.7.13

python /scratch/groups/khosla/Orphan_PKS/scripts/counting_domains_sk.py -cluster_file1 /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundant.alldbs.txt -cluster_file2 /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt -inputf /scratch/groups/khosla/Orphan_PKS/results/antismash_results/alldbs -outputf /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs
```
Results go to clusterTableNonRedundantNonSimilar.alldbs.prot_len.correct_domains.txt

**15. Creating the distance matrix**
This script reads all similarity scores and makes a distance matrix.
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/pairsToMatrix_blastp.pl /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.correct_domains.txt> /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/dataMatrix_alldbs_blastp.txt
```

***16. Creating the catalog of NRNSNSS***
NOTE: In cases where a cluster that was deemed redundant to another main one, had identical clusters (=same sequence or same species+architecture), then these clusters were also added as redundant to the main cluster, into the last column.
```
python /scratch/groups/khosla/Orphan_PKS/scripts/makeClusterTableSummaryNonSequenceSimilar_sk.py -cluster_tableNRNS /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.correct_domains.txt -nonredundant_file /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs.nonredundant_chained_new.cutoff0.9.txt -redundant_set_file /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs.redundant_chained_new.cutoff0.9.txt -output_f /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs
```
The result goes to:
clusterTableNonRedundantNonSimilarNonSequenceSimilar.alldbs.prot_len.correct_domains.txt

#########################################
#########################################
### Known Cluster Finding:
We want to extract PKS-containing clusters from MIBig, which wee consider as "known" clusters. <br />
We then look through our data and try to find whch clusters correspond to these "known" ones. <br />
We also want to look at cluster annotations, where some of them are moted as "known" clusters. <br />

### Finding known PKS clusters in MIBig
First, we look at MIBig 3.1 json files to extract all clusters: not only PKS, but also  NRPS, RiPPs etc. <br />
We will match relevant clusters to our list of PKSs later.

Download MiBIG Product Files from here: https://mibig.secondarymetabolites.org/download <br />
Downloaded 3.1 JASON and GBK files to local computer <br />
```
module load python/2.7.13
module load biology py-biopython/1.70_py27
module load math py-numpy
python extracting_all_from_json_sk.py
```
The output goes to: all_clusters_mibig.txt

### Matching known clusters from MIBig to our list of clusters

```
python finding_all_correspondencies_from_MIBig_sk.py 
```
The result goes to: <br />
known_all_clusters_from_MIBig.txt <br />
There is a total of 339 clusters matched. <br />

### Looking for clusters annotated as known 
For some PKS clusters, their description mentions if they are known. We want to list known clusters matched to to clusters from the ClusterTable (before finding redundants)

```
python finding_PKS_correspondencies_from_clusterTable_sk.py
```
The result goes to: <br />
known_clusters_from_clusterTable.txt <br />
Some of them have to be manually modified, because it's difficult to correctly parse all product names. <br />

The result goes to: <br />
known_clusters_from_clusterTable_edited.txt <br />
There is a total of 426 clusters identified this way. Some of them overlap with the list extracted from MIBig. <br />

### Assembling the list of knowm clusters from both sources

```
python assembling_correspondencies.py
```
The results go to: <br />
known_clusters_MIBig_and_clusterTable.txt <br />
There is a total of 508 known clusters. <br />
Then I manually edit for typos etc. <br />
The result goes to: # known_clusters_MIBig_and_clusterTable_edited.txt <br />

### Matching this list to clusters that appear in the dendrogram (NRNS) 

```
module load python/2.7.13
python matching_correspondencies_to_all_leaves_sk.py 
```
The result goes to: # known_clusters_MIBig_and_clusterTable_matched_to_all_leaves.txt <br />
There is a total of 492 known clusters matched to clusters from the clusterTableNRNS. <br />

### Matching this list to "main" clusters in the dendrogram (NRNSNSS)

```
module load python/2.7.13
python matching_correspondencies_to_main_leaves_sk.py 
```
The result goes to: # known_clusters_MIBig_and_clusterTable_matched_to_main.txt <br />
There is a total of 437 known clusters matched to clusters from the clusterTableNRNS. <br />

### Making a timeline

### ClusterTables for different years 
```
module load python/2.7.13
module load math py-numpy
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ClusterTablesTimeline
python /scratch/groups/khosla/Orphan_PKS/scripts/making_clusterfiles_by_year_sk.py
```

### DistanceMatrices for different years 

```
sbatch creating_distance_matrices_example.sh
```

### ParsedScores for different years

```
module load python/2.7.13
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores
sbatch generating_parsedscores.sh
```
Runs generating_parsedscore_by_year_sk.py

### Looking for redundancies in scores and chaining for different years 

```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores
sbatch chaining_1.sh  
```
runs chaining_blastp.pl

### Making graphs for the number of non-redundant clusters, and the percentage of clusters with redundancies 

```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores
sbatch creating_timeline.sh   (runs creating_timeline_sk.py)
```
The results go to: <br />
timeline_unique_red.csv (numbers) <br />
Timeline_unique_redundant_clusters.png (graph) (plots the number of unique clusters per year)<br />

Opened the csv file in excel and used prism to plot the rediscovery rate i.e. number of redundant clusters divided by total number of clusters

### Plotting the dendrogram 
This python script takes the distance matrix file, the list of redundant and non-redundant clusters, and performs clustering only for distinct entries.  As in the 2013 paper,we are using the McQuitty method (also called "weighted" in scipy) for clustering. Clusters with redundancies are marked "seen more than once". This script shows names of known clusters and highlights them in red on the dendrogram of either all clusters, or distinc clusters.  The list of known distinct clusters is: known_clusters_MIBig_and_clusterTable_matched_to_main_edited.txt 


```
module load python/2.7.13
module load biology py-biopython/1.70_py27
module load math py-numpy
ml math py-scipy
ml math py-pandas
ml math py-sympy
ml viz py-matplotlib
module load python/2.7.13
python plotting_dendrograms_sk.py -show_all No -show_known Yes -only_show_known No -show_accession Yes -show_species Yes -show_date No -show_type Yes -show_domains No -color_known Yes
```

### How similar are orphan clusters from 2022 to known clusters?

```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/Similarity
sbatch similarity_to_known_clusters.sh (runs similarity_to_known_clusters_sk.py)
```

### How similar are distinct orphan clusters from 2022 to distinct known clusters?

```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/Similarity
sbatch similarity_to_known_clusters_distinct.sh (runs Similarity_to_known_clusters_distinct_sk.py)
```

### How similar are orphan clusters to known clusters, along the years?
```
cd /scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/Similarity
sbatch similarity_to_known_clusters_distinct_timeline.sh (runs similarity_to_known_clusters_distinct_timeline_sk.py)
```
