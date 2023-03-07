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

with open("known_clusters_MIBig_and_clusterTable_reedited.txt",'r') as known_cluster_file:
	known_clusters=known_cluster_file.readlines()

known_assembly_line_PKS_dict={}
check_appeared={}
for line in known_clusters:
	line_info=line.split("\t")
	cluster_id=line_info[0].strip(":")
	product=line_info[1].strip("\n")
	known_assembly_line_PKS_dict[cluster_id]=product
	check_appeared[cluster_id]=0

# print "LC330869 1",known_assembly_line_PKS_dict["LC330869 1"]
with open("known_clusters_MIBig_and_clusterTable_matched_to_all_leaves.txt",'w') as output_file:
	with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.correct_domains.txt",'r') as final_cluster_file:
		final_cluster_list=final_cluster_file.readlines()
        # Reading the final cluster list and splitting each line into a list,
        # where each element corresponds to a cluster with its description.
        # Then, splitting each element into another list,
        # where each element corresponds to a descriptor (accession, cluster name, strain etc)
        for n in range(1,len(final_cluster_list)):
            line=final_cluster_list[n].strip("\n")
            line=line.split("\t")   
            cluster_ids=[]

            # Now we can look at the descriptors for the cluster that serves as the main accession
            main_accession=line[0]
            main_cluster_n=str(line[4][9:-1])
            cluster_ids.append(main_accession+" "+main_cluster_n)

            # Now we can look at the descriptors for the clusters that are considered redundant to the main accession

            addiitonal_clusters1=line[11].split("; ")
            for m in range(0,len(addiitonal_clusters1)):
                redund_accession_data=addiitonal_clusters1[m].split(" ")
                # There are cases when the current cluster data is too short to extract data from:
                # it happens when ";" appears in the desciption, and part of the description is considered to describe another cluster, which it does not.
                # We do nothing about these cases, just care about the correct ones.
                if len(redund_accession_data)>1:
                    if redund_accession_data[1][0:5]=="Date=":
                        accession=redund_accession_data[0]
                        if redund_accession_data[2]=="Cluster":
                            cluster_n=redund_accession_data[3][1:-1]
                            if cluster_n[0:3]!="r1c":
                                print "Not the right format",cluster_n
                            cluster_ids.append(accession+" "+cluster_n)
                        elif redund_accession_data[2][0:7]=="cluster":
                            cluster_n=redund_accession_data[2][8:-1]
                            if cluster_n[0:3]!="r1c":
                            	print "Not the right format",cluster_n
                            cluster_ids.append(accession+" "+cluster_n)
                        else:
                            print "No cluster number found for ",accession,clusternum

            addiitonal_clusters2=line[12].split("; ")
            for m in range(0,len(addiitonal_clusters2)):
                redund_accession_data=addiitonal_clusters2[m].split(" ")
                if len(redund_accession_data)>1:
                    # There are cases when the current cluster data is too short to extract data from:
                    # it happens when ";" appears in the desciption, and part of the description is considered to describe another cluster, which it does not.
                    # We do nothing about these cases, just care about the correct ones.
                    if redund_accession_data[1][0:5]=="Date=":
                        accession=redund_accession_data[0]
                        if redund_accession_data[2]=="Cluster":
                            cluster_n=redund_accession_data[3][1:-1]
                            if cluster_n[0:3]!="r1c":
                                print "Not the right format",cluster_n
                            cluster_ids.append(accession+" "+cluster_n)
                        elif redund_accession_data[2][0:7]=="cluster":
                            cluster_n=redund_accession_data[2][8:-1]
                            if cluster_n[0:3]!="r1c":
                                print "Not the right format",cluster_n

                            cluster_ids.append(accession+" "+cluster_n)
                        else:
                            print "No cluster number found for ",accession,clusternum

            # addiitonal_clusters3=line[13].split("; ")
            # for m in range(0,len(addiitonal_clusters3)):
            #     redund_accession_data=addiitonal_clusters3[m].split(" ")
            #     if len(redund_accession_data)>1:
            #         # There are cases when the current cluster data is too short to extract data from:
            #         # it happens when ";" appears in the desciption, and part of the description is considered to describe another cluster, which it does not.
            #         # We do nothing about these cases, just care about the correct ones.
            #         if redund_accession_data[1][0:5]=="Date=":
            #             accession=redund_accession_data[0]
            #             if redund_accession_data[2]=="Cluster":
            #                 cluster_n=redund_accession_data[3]
            #                 cluster_ids.append(accession+" "+cluster_n)
            #             elif redund_accession_data[2][0:7]=="cluster":
            #                 cluster_n=redund_accession_data[2][7:]
            #                 cluster_ids.append(accession+" "+cluster_n)
            #             else:
            #                 print "No cluster number found for ",accession,clusternum
            
            if main_accession+" "+main_cluster_n=="LC330869 1":
                print "current ids:",cluster_ids
           
            #print "current ids:",cluster_ids 
            known_products=[]
            # Going through each cluster ID in this line
            for curr_cluster in cluster_ids:

                # checking if there is a known assembly-line PKS or other cluster that overlaps with this location:
                if (curr_cluster in known_assembly_line_PKS_dict) and (known_assembly_line_PKS_dict[curr_cluster] not in known_products):
                	known_products.append(known_assembly_line_PKS_dict[curr_cluster])
		if curr_cluster in known_assembly_line_PKS_dict and curr_cluster in check_appeared:
			del check_appeared[curr_cluster]
            if len(known_products)!=0:
	            output_file.write(main_accession+" "+main_cluster_n+"\t"+", ".join(known_products)+"\n")

print "Printing accessions that didn't appear",check_appeared
