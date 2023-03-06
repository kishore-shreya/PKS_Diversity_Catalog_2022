# PKS_Diversity_Catalog_2022

########################################################################################
#########################################################################################
 
                                 DISTINCT NON-REDUNDANT PKS CATALOG (2022)

########################################################################################
#########################################################################################

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
Archae: 2157 <br />
Bacteria: 2 <br />
Eukaryota: 2759 <br />

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

