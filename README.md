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

Repeat for the remaining 7 databases

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
Repeat for the remaining 7 databases

######################################################################################################
# *****
### Filter genomes with at least 3 KSs
```
perl /scratch/groups/khosla/Orphan_PKS/scripts/filterGenomes.pl KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart 3 > KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3 
```
Repeat for the remaining 7 databases and SRA WGS

### How many entries meet the criteria?

cut -f 1 KSSignature_clusters.patnt.max.20000.min.3000.bp.apart.min3 |sort -u |wc <br />
383  383 4191 <br />

Repeat for the remaining 7 databases and SRA WGS

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

Repeat for the remaining 7 databases and SRA WGS

### Separate Genbank File into Individual Files
```
ml perl
cpanm install Bio::SeqIO
perl /scratch/groups/khosla/Orphan_PKS/scripts/separateSequenceFiles.pl /scratch/groups/khosla/Orphan_PKS/KS_Cluster_results/KSSignature_clusters.env_nt.max.20000.min.3000.bp.apart.min3.gbwithparts genbank genbank.env_nt.min3.newset1
```
creates new folder called "genbank.env_nt.min3" and puts all individual files in there <br />
Repeat for the remaining 7 databases and SRA WGS

### Convert Genbanks into FASTA

Need to run this in the genbank.env_nt.min3.newset1 directory
```
ml perl
cpanm install Bio::SeqIO
for i in *.gb; do
perl /scratch/groups/khosla/Orphan_PKS/scripts/gb2fasta.pl /scratch/groups/khosla/Orphan_PKS/KS_Cluster_results/genbank.env_nt.min3.newset1/$i /scratch/groups/khosla/Orphan_PKS/KS_Cluster_results/env_nt_fasta/${i%%.*}.fasta
done
```
Repeat for the remaining 7 databases and SRA WGS

##########################################################################################################################
# *****

### Install AntiSMASH Locally on Sherlock




