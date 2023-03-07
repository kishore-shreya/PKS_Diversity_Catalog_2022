#!/usr/local/bin/python

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

#Shreya Kishore edits to make it work with antiSMASH version 6 output files instead of antiSMASH version 4 that was previously used



# This script  that:
# - for the clusters listed in that clusterTable
# - retrieves their genbank files
# - extracts the protein sequences of biosynthetic enzymes from these clusters
# - and saves them into fasta files.

# It takes as arguments:
# - the clusterTable file
# - the range of clusters to process (e.g. 0 1000; or 1000 2000 )
# - the location of input genbank files, without the subfolders (e.g. antismash_from_fasta)
# - the location of output fasta files

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
import os


parser = argparse.ArgumentParser(description='Get the range of cluster IDs to parse')

parser.add_argument('-cluster_file','--clusterfile')
parser.add_argument('-startID','--start_ID')
parser.add_argument('-endID','--end_ID')
parser.add_argument('-inputf','--input_folder')
parser.add_argument('-outputf','--output_folder')

args=parser.parse_args()

cluster_file=args.clusterfile
cluster_ID_start=int(args.start_ID)
cluster_ID_end=args.end_ID
my_in_path=args.input_folder
my_out_path=args.output_folder

with open(cluster_file,'r') as cluster_list_file:
    cluster_list=cluster_list_file.readlines()

cluster_ids=[]
for item in cluster_list:
    info=item.split("\t")
    accession=info[0]
    cluster_n=info[4][8:]
    cluster_ids.append((accession,cluster_n))

cluster_ids=cluster_ids[1:]
if cluster_ID_end=="end":
    cluster_ID_end=len(cluster_ids)
else:
    cluster_ID_end=int(cluster_ID_end)
cluster_ids_curr=cluster_ids[cluster_ID_start:cluster_ID_end]


with open("List_of_missing_clusterfiles.txt",'w') as missing_cluster_file:
    for item in cluster_ids_curr:
        accession=item[0]
        cluster_num=item[1]
	cluster_length=len(cluster_num)
	#print "length of this cluster =", cluster_length
        if len(cluster_num)==6:
             cluster_n="region00"+cluster_num[4]
        elif len(cluster_num)==7:
	     cluster_n="region0"+cluster_num[4:6]
        else:
             print "Error: cluster number of unexpected size, "+accession+" "+cluster_num
        subfoldername="/"+accession
	
	filename="/"+accession+"."+cluster_n
	filename2="/"+accession+".1."+cluster_n

	path = my_in_path+subfoldername+filename+".gbk"
	path2 = my_in_path+subfoldername+filename2+".gbk"
	if (os.path.exists(path)):
     	    gbk_file = path
	elif (os.path.exists(path2)):
	    gbk_file = path2

        try:
            open(gbk_file,'r')
        except:
            missing_cluster_file.write(gbk_file+"\n")
        else:
            print "opened genbank file"
            record = SeqIO.read(gbk_file, "gb")
            feature_list=record.features
            feature_n=len(feature_list)
            protein_seq=Seq("",IUPAC.protein)
        
            with open(my_out_path+filename+".fasta", "w") as output_file:
		print my_out_path+filename+".fasta"

                # going over CDS features
                for i1 in range(0,feature_n):
                    feature1=feature_list[i1]
                    if feature1.type=="CDS":

                        # going over aSDomain features within that CDS
                        for i2 in range(i1+1,feature_n):
                            feature2=feature_list[i2]
                            if feature2.location.start in feature1 and (feature2.type=="aSDomain"):
                                # print feature1.qualifiers["translation"][0]
                                protein_seq+=feature1.qualifiers["translation"][0]
                                break;
                
                record_id=record.id+"."+cluster_num
                
                # save this record into a fasta file
                protein_record=SeqRecord(protein_seq, id=record_id, description=record.description)
                SeqIO.write(protein_record, output_file, "fasta")
