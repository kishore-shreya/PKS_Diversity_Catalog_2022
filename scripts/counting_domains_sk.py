#!/usr/local/bin/python

#shreya kishore edits to make it work with new antismash 6 file structure

# Copyright (C) 2019 Aleksandra Nivina

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# To see the GNU General Public License, Plesase see 
# <http://www.gnu.org/licenses/>.

# This script was modified relative to the initial one.
# Before, there were two scripts: one that found redundancies, and another that chained them.
# Now, there's only one script that finds redundancies and chains them, in a smarter way.
# The resulting files are two lists:
# - list of non-redundant clusters, with 1 cluster/line
# - list of redundant clusters, where the first cluster/line is the "main" record, and the others are redundant
# These lists are stored in the same place as the input file with similarity scores.

# This script:
# - for the clusters listed in that clusterTable
# - retrieves their genbank files
# - extracts the number of each type of enzymatic domain
# - records them into new cluster tables


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse


parser = argparse.ArgumentParser(description='Get the range of cluster IDs to parse')

parser.add_argument('-cluster_file1','--clusterfile1')
parser.add_argument('-cluster_file2','--clusterfile2')
parser.add_argument('-inputf','--input_folder')
parser.add_argument('-outputf','--output_folder')

args=parser.parse_args()

cluster_file1=args.clusterfile1
cluster_file2=args.clusterfile2
my_in_path=args.input_folder
my_out_path=args.output_folder


############## Dealing with the clusterTable NR ##############

# clusterTable NR
with open(cluster_file1,'r') as cluster_list_file1:
    cluster_list1=cluster_list_file1.readlines()

# Changing the table header
header1=cluster_list1[0].split("\t")
del header1[10:18]
header1[9]="List of domains"

# Changing the rest of the table, and saving the domain number description into a dictionary
domain_numbers={}
with open(my_out_path+"/clusterTableNonRedundant.alldbs.correct_domains.txt",'w') as output_file1:
    output_file1.write("\t".join(header1))

    for i in range(1,len(cluster_list1)):
        line=cluster_list1[i]
        info= line.split("\t")
        accession=info[0]
        cluster_num=info[4][9:-1]
	#print "cluster_num is", cluster_num
	#print "length is", len(cluster_num)
        if len(cluster_num)==4:
            cluster_n="00"+cluster_num[3]
        elif len(cluster_num)==5:
            cluster_n="0"+cluster_num[3:]
        elif len(cluster_num)==6:
            cluster_n=cluster_num[3:]
        subfoldername="/"+accession+"/"
	#print "cluster_n is", cluster_n
        filename=accession+".region"+cluster_n

        n_CAL=0
        n_KS=0
        n_AT=0
        n_GNAT=0
        n_tATd=0
        n_KR=0
        n_DH=0
        n_ER=0
        n_Hal=0
        n_ECH=0
        n_PS=0
        n_ACP=0
        n_ACPb=0
        n_CYC=0

        n_C=0
        n_A=0
        n_E=0
        n_F=0
        n_H=0
        n_PCP=0

        n_AmT=0
        n_MT=0
        n_FkbH=0
        n_TE=0
        n_TD=0
        n_NADb=0

	try:
		record = SeqIO.read(my_in_path+subfoldername+filename+".gbk", "gb")

	except IOError:

		filename1=accession+".1.region"+cluster_n
		record = SeqIO.read(my_in_path+subfoldername+filename1+".gbk", "gb")

        for feature in record.features:
            if feature.type=="aSDomain":



                # PKS domains
                if feature.qualifiers.get("aSDomain")[0]=="CAL_domain":
                    n_CAL+=1

                if feature.qualifiers.get("aSDomain")[0]=="PKS_KS":
                    n_KS+=1
  
                elif feature.qualifiers.get("aSDomain")[0]=="PKS_AT":
                    n_AT+=1

                elif feature.qualifiers.get("aSDomain")[0]=="GNAT" :
                    n_GNAT+=1

                elif feature.qualifiers.get("aSDomain")[0]=="Trans-AT_docking":
                    n_tATd+=1

                elif feature.qualifiers.get("aSDomain")[0]=="PKS_KR":
                    n_KR+=1

                elif feature.qualifiers.get("aSDomain")[0]=="PKS_DH" or feature.qualifiers.get("aSDomain")[0]=="PKS_DH2" or feature.qualifiers.get("aSDomain")[0]=="PKS_DHt":
                    n_DH+=1


                elif feature.qualifiers.get("aSDomain")[0]=="PKS_ER":
                    n_ER+=1

                elif feature.qualifiers.get("aSDomain")[0]=="Hal":
                    n_Hal+=1

                elif feature.qualifiers.get("aSDomain")[0]=="ECH":
                    n_ECH+=1

                elif feature.qualifiers.get("aSDomain")[0]=="PS":
                    n_PS+=1

                elif feature.qualifiers.get("aSDomain")[0]=="ACP":
                    n_ACP+=1

                elif feature.qualifiers.get("aSDomain")[0]=="ACP_beta" :
                    n_ACPb+=1

                elif feature.qualifiers.get("aSDomain")[0]=="Polyketide_cyc" :
                    n_CYC+=1



                #NRPS domains
                elif feature.qualifiers.get("aSDomain")[0]=="Condensation" or feature.qualifiers.get("aSDomain")[0]=="Condensation_DCL" or feature.qualifiers.get("aSDomain")[0]=="Condensation_Dual" or feature.qualifiers.get("aSDomain")[0]=="Condensation_LCL" or feature.qualifiers.get("aSDomain")[0]=="Condensation_Starter":
                    n_C+=1

                elif feature.qualifiers.get("aSDomain")[0]=="AMP-binding":
                    n_A+=1

                elif feature.qualifiers.get("aSDomain")[0]=="Epimerization*":
                    n_E+=1

                elif feature.qualifiers.get("aSDomain")[0]=="F":
                    n_F+=1

                elif feature.qualifiers.get("aSDomain")[0]=="Heterocyclization":
                    n_H+=1  

                elif feature.qualifiers.get("aSDomain")[0]=="PCP":
                    n_PCP+=1



                # PKS/NRPS domains
                elif feature.qualifiers.get("aSDomain")[0]=="Aminotran_1_2" or feature.qualifiers.get("aSDomain")[0]=="Aminotran_3" or feature.qualifiers.get("aSDomain")[0]=="Aminotran_4" or feature.qualifiers.get("aSDomain")[0]=="Aminotran_5":
                    n_AmT+=1

                elif feature.qualifiers.get("aSDomain")[0]=="MT" or feature.qualifiers.get("aSDomain")[0]=="cMT" or feature.qualifiers.get("aSDomain")[0]=="nMT" or feature.qualifiers.get("aSDomain")[0]=="oMT":
                    n_MT+=1

                elif feature.qualifiers.get("aSDomain")[0]=="FkbH":
                    n_FkbH+=1

                elif feature.qualifiers.get("aSDomain")[0]=="Thioesterase":
                    n_TE+=1

                elif feature.qualifiers.get("aSDomain")[0]=="TD":
                    n_TD+=1

                elif feature.qualifiers.get("aSDomain")[0]=="NAD_binding_4":
                    n_NADb+=1
                    
        domain_number_list=[n_CAL,n_KS,n_AT,n_GNAT,n_tATd,n_KR,n_DH,n_ER,n_Hal,n_ECH,n_PS,n_ACP,n_ACPb,n_CYC,n_C,n_A,n_E,n_F,n_H,n_PCP,n_AmT,n_MT,n_FkbH,n_TE,n_TD,n_NADb]
        domain_name_list=["CAL","KS","AT","GNAT","tATd","KR","DH","ER","Hal","ECH","PS","ACP","ACPb","CYC","C","A","E","F","H","PCP","AmT","MT","FkbH","TE","TD","NADb"]

        # deleting columns that corresponded to domain numbers, and creating one column with text description of all detected domains
        del info[10:18]
        curr_domains=""
        # print "Found some domains:",n_KS,n_AT
        for n in range(0,len(domain_number_list)):
            # print domain_name_list[n],":",domain_number_list[n]
            if domain_number_list[n]>0:
                curr_domains+=" "+str(domain_number_list[n])+domain_name_list[n]
        info[9]=curr_domains


        domain_numbers[(accession,cluster_num)]=curr_domains

        new_line="\t".join(info)
        output_file1.write(new_line)



############## Dealing with the clusterTable NRNS ##############


# clusterTable NRNS
with open(cluster_file2,'r') as cluster_list_file2:
    cluster_list2=cluster_list_file2.readlines()

# Changing the table header
header2=cluster_list2[0].split("\t")
del header2[10:18]
header2[9]="List of domains"

# Changing the rest of the table, using the information saved in the dictionary
with open(my_out_path+"/clusterTableNonRedundantNonSimilar.alldbs.prot_len.correct_domains.txt",'w') as output_file2:
    output_file2.write("\t".join(header2))
    for i in range(1,len(cluster_list2)):
        line=cluster_list2[i]
        info= line.split("\t")
        accession=info[0]
        cluster_num=info[4][9:-1]

        curr_domains=domain_numbers[(accession,cluster_num)]

        del info[10:18]
        info[9]=curr_domains
        new_line="\t".join(info)
        output_file2.write(new_line)



# Number of different domains detected
# (description can be found at https://docs.antismash.secondarymetabolites.org/modules/nrps_pks_domains/):

# Aminotran_5 717
# GNAT 499
# PKS_DH 25535
# Hal 13
# Epimerization 2803
# F 26
# PKS_Docking_Nterm 9567
# Thioesterase 10141
# B 84
# ACP_beta 4134
# PS 137
# TD 768
# AMP-binding 22126
# Trans-AT_docking 18652
# Condensation 17402
# NAD_binding_4 1227
# FkbH 1220
# CAL_domain 3070
# PKS_ER 10161
# NRPS-COM_Cterm 269
# PKS_ACP 0
# Heterocyclization 2028
# Polyketide_cyc 676
# PKS_Docking_Cterm 11876
# ACPS 1688
# Aminotran_3 2033
# PKS_AT 51451
# ECH 2832
# MT 7568
# ACP 75027
# PKS_KR 60582
# PKS_KS 79263
# NRPS-COM_Nterm 882
# PKS_DH2 13567
# Aminotran_4 179
# PKS_DHt 2825
# PCP 24736
# Aminotran_1_2 2297      
