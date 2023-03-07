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

import numpy as np
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# reading the list of known assembly-line PKS clusters and saving their IDs and coordinates
# key = accession; value = [(start,end,MIBig_accession,Product),(start,end,MIBig_accession,Product),...]
# with open("/scratch/groups/khosla/Orphan_PKS/10MarNonOrphans/assembly_line_pks_mibig.txt",'r') as known_assembly_line_file:
with open("./all_clusters_mibig.txt",'r') as known_assembly_line_file:
    known_assembly_line_PKS_list=known_assembly_line_file.readlines()
known_assembly_line_PKS_dict={}
for line in known_assembly_line_PKS_list:
    line=line.split(",")
    # accession=line[0]
    MIBig=line[1]
    product=line[2].strip("\n")
    if product[0].islower():
        product=product.capitalize()
    if MIBig!="BGC0001129":
    # for record in SeqIO.parse("/scratch/groups/khosla/Orphan_PKS/10MarNonOrphans/mibig_gbk_1.4/"+MIBig+".gbk", "gb"):
        for record in SeqIO.parse("./mibig_gbk_1.4/"+MIBig+".gbk", "gb"):

            comment=record.annotations["comment"]
            if comment.find("GenBank ID")!=-1:
                accession_start=comment.find("GenBank ID")+11
                if comment.find(".",accession_start)!=-1:
                    accession=comment[accession_start:comment.find(".",accession_start)]
                else:
                    accession=comment[accession_start:].strip("\n").strip(" ")
                # if the cluster is a subset of that accession, find locations here
                if comment.find("region between ")!=-1:
                    cluster_location_start_pos=comment.find("region between ")+15
                    cluster_location_separator_pos=comment.find("-",cluster_location_start_pos)
                    cluster_location_end_pos=comment.find("nt",cluster_location_start_pos)
                    cluster_start=int(comment[cluster_location_start_pos:cluster_location_separator_pos])
                    cluster_end=int(comment[cluster_location_separator_pos+1:cluster_location_end_pos])
                    if accession not in known_assembly_line_PKS_dict:
                        known_assembly_line_PKS_dict[accession]=[(cluster_start,cluster_end,MIBig,product)]
                        break
                    else:
                        known_assembly_line_PKS_dict[accession].append((cluster_start,cluster_end,MIBig,product))
                        break

                # otherwise, locations is the entire accession
                else:
                    for feature in record.features:
                        if feature.type=="source":
                            cluster_start=int(feature.location.start)
                            cluster_end=int(feature.location.end)
                            if accession not in known_assembly_line_PKS_dict:
                                known_assembly_line_PKS_dict[accession]=[(cluster_start,cluster_end,MIBig,product)]
                                break
                            else:
                                known_assembly_line_PKS_dict[accession].append((cluster_start,cluster_end,MIBig,product))
                                break



    # reading the cluster file with non-identical records, and making a dictionary of cluster cordinates for each cluster ID.
    # key = (accession,cluster_number); value = (start,end)
    # with open("/scratch/groups/khosla/Orphan_PKS/10MarNonOrphans/clusterTableNonRedundant.6dbs.txt",'r') as nonidential_cluster_file:
    with open("known_all_clusters_from_MIBig.txt",'w') as output_file1: 

        with open("./clusterTableNonRedundant.6dbs.correct_domains.txt",'r') as nonidential_cluster_file:
            nonidentical_cluster_list=nonidential_cluster_file.readlines()
        nonidentical_cluster_list=nonidentical_cluster_list[1:]
        # cluster_locations={}
        for item in nonidentical_cluster_list:
            item_data=item.split("\t")
            accession=item_data[0]
            # print item_data[4]
            cluster_n=int(item_data[4][7:])
            coordinates=item_data[6]
            clusterstart=int(coordinates[:coordinates.find("-")])
            clusterend=int(coordinates[coordinates.find("-")+1:])
            # cluster_locations[(accession,cluster_n)]=(clusterstart,clusterend)
            if accession in known_assembly_line_PKS_dict :
                known_clusters=known_assembly_line_PKS_dict[accession]
                for known_cluster in known_clusters:
                # print known_cluster
                    known_cluster_start=known_cluster[0]
                    known_cluster_end=known_cluster[1]
                    known_cluster_MIBig=known_cluster[2]
                    known_cluster_product=known_cluster[3]
                    if (clusterstart<=known_cluster_start and clusterend>=known_cluster_end)\
                    or (clusterstart>=known_cluster_start and clusterend<=known_cluster_end)\
                    or (clusterstart<known_cluster_end and clusterend>=known_cluster_end)\
                    or (clusterstart<=known_cluster_start and clusterend>known_cluster_start):
                        # known_matches.append((known_cluster_MIBig,known_cluster_product))
                        # print "matched",curr_accession,curr_start,curr_end, "with", known_cluster
                        output_file1.write(accession+" "+str(cluster_n)+":\t"+known_cluster_product+"\n")
