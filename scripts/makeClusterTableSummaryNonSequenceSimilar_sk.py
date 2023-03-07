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


# This script takes the list of clusters that were compared by plastp (cluster_tableNRNS),
# the list of redundant sets and the list of nonredundant clusters,
# and creates a new clusterTableNRNSNSS that has only the main leaves listed.
# Its last 3 columns correspond to redundant, similar and sequence similar clusters.

# python makeClusterTableSummaryNonSequenceSimilar.py \
# -cluster_tableNRNS ../11_counting_domains/clusterTableNonRedundantNonSimilar.6dbs.prot_len.correct_domains.txt  \
# -nonredundant_file ../12_plotting_dendrograms/parsed_scores.all_final.nonredundant_chained_new.cutoff0.9.txt \
#  -redundant_set_file ../12_plotting_dendrograms/parsed_scores.all_final.redundant_chained_new.cutoff0.9.txt \
#  -output_f ./

import argparse

parser = argparse.ArgumentParser(description='Get input files')

parser.add_argument('-cluster_tableNRNS','--clustertableNRNS')
parser.add_argument('-nonredundant_file','--nonredundantfile')
parser.add_argument('-redundant_set_file','--redundantsetfile')
parser.add_argument('-output_f','--outputfolder')

args=parser.parse_args()

cluster_tableNRNS=args.clustertableNRNS
nonredundant_file=args.nonredundantfile
redundant_set_file=args.redundantsetfile
output_folder=args.outputfolder


# Reading the NRNS cluster table file and creating a short description for all clusters,
# and a list of redundant cluster IDs and descriptions.
# Many cluster IDs in the previosuly listed redundant lists are repeated. Fixing that here.
short_cluster_descr={}
previous_cluster_redundants1={}
previous_cluster_redundants2={}
redundant_cluster_short_descriptions={}
# redundant_cluster_short_descriptions_aa={}
with open(cluster_tableNRNS,'r') as cluster_list_file:
    cluster_list=cluster_list_file.readlines()
for i in range(1,len(cluster_list)):
    line=cluster_list[i]
    info= line.strip("\n").split("\t")
    accession=info[0]
    cluster_num=info[4][9:-1].strip(" ")
    cluster_id=accession+" "+cluster_num
    # print cluster_id
    date=info[2]
    length=info[7]
    description=info[1]
    redundants1=info[11].split("; ")
    redundants2=info[12].split("; ")
    short_descr=accession+" Date="+date+" cluster"+cluster_num+" "+length+"aa "+description.strip("\n")
    # redundant_cluster_short_descriptions_aa[cluster_id]=short_descr
    short_cluster_descr[cluster_id]=short_descr
    # redundants=redundants1+redundants2
    empty_item=""
    while empty_item in redundants1:
		redundants1.remove(empty_item)

    redundant_id_list1=[]
    for redund in redundants1:
		redund_info=redund.strip("\n").split(" ")
		if len(redund_info)>2:
			redund_accession=redund_info[0]
			if redund_info[2]=="Cluster" or redund_info[2]=="cluster":
				redund_clusternum=redund_info[3]
			else:
				redund_clusternum=redund_info[2][7:]
			redund_cluster_id=redund_accession+" "+redund_clusternum

			# recording a non-previously seed cluster ID among redundants, and recording a short description for it
			if (redund_cluster_id!=cluster_id) and (redund_cluster_id not in redundant_id_list1):
				# print "main:",cluster_id,"redundant:",redund_cluster_id
				redundant_id_list1.append(redund_cluster_id)
				redundant_cluster_short_descriptions[redund_cluster_id]=redund.strip("\n")
				
		else:
			print "short entry in the first loop:",redund_info

    previous_cluster_redundants1[cluster_id]=redundant_id_list1

    while empty_item in redundants2:
		redundants2.remove(empty_item)
    redundant_id_list2=[]
    for redund in redundants2:
		redund_info=redund.strip("\n").split(" ")
		if len(redund_info)>2:
			redund_accession=redund_info[0]
			if redund_info[2]=="Cluster" or redund_info[2]=="cluster":
				redund_clusternum=redund_info[3]
			else:
				redund_clusternum=redund_info[2][7:]
			redund_cluster_id=redund_accession+" "+redund_clusternum

			# recording a non-previously seed cluster ID among redundants, and recording a short description for it
			# recording a short description for each cluster previously recorded as redundant
			if (redund_cluster_id!=cluster_id) and (redund_cluster_id not in redundant_id_list1) and (redund_cluster_id not in redundant_id_list2):
				# print "main:",cluster_id,"redundant:",redund_cluster_id
				redundant_id_list2.append(redund_cluster_id)
				redundant_cluster_short_descriptions[redund_cluster_id]=redund.strip("\n")
		else:
			print "short entry in the first loop:",redund_info

    previous_cluster_redundants2[cluster_id]=redundant_id_list2

    # if cluster_id=="PHNC01000011 1":
    # 	print cluster_id
    # 	print "Redundant:",redundant_id_list1
    # 	print "Similar:",redundant_id_list2



# Reading the list of non-redundant clusters and creating a dictionary
non_redundant_clusters={}
with open(nonredundant_file,'r') as non_redundant_file:
	nonredundant_list=non_redundant_file.readlines()
for line in nonredundant_list:
	nonredund_cluster=line.strip("\n")
	non_redundant_clusters[nonredund_cluster]=1

# Reading the list of redundant clusters and creating 2 dictionaries that match main to redundants and each redundant to mains (cluster IDs)
main_to_redundant_clusters={}
redundant_to_main_clusters={}
with open(redundant_set_file,'r') as redundant_set_file:
	redundant_set_list=redundant_set_file.readlines()
for line in redundant_set_list:
	redund_clusters=line.strip("\n").split("\t")
	main_cluster=redund_clusters[0]
	redund_clusters=redund_clusters[1:]

	# deleting empty items fromt he list of redundant clusters
	empty_item=""
	while empty_item in redund_clusters:
		redund_clusters.remove(empty_item)

	# recording into dictionary that matches main entries with all redundant clusters (cluster IDs)
	main_to_redundant_clusters[main_cluster]=redund_clusters

	# if main_cluster=="PHNC01000011 1":
	# 	print "\nfound a list of redundant clusters:",redund_clusters

	# recording into dictionary that matches each redundant cluster with all main clusters (cluster IDs)
	for redund_cluster in redund_clusters:
		redund_cluster=redund_cluster.strip("\n")
		if redund_cluster in redundant_to_main_clusters:
			redundant_to_main_clusters[redund_cluster].append(main_cluster)
			# if main_cluster=="PHNC01000011 1":
			# 	print redund_cluster,"already accounted for"
		else:
			redundant_to_main_clusters[redund_cluster]=[main_cluster]
			# if main_cluster=="PHNC01000011 1":
			# 	print redund_cluster,"not yet accounted for"

# print "For Redundant AM778535 1, the redundant clusters are:",redundant_to_main_clusters["AM778535 1"]
# print "For Main PHNC01000011 1, the redundant clusters are:",main_to_redundant_clusters["PHNC01000011 1"]

# Recording the new clusterTable file

# Generating additional descriptions

# This dictionary will keep Cluster IDs of all additional redundant clusters
main_to_additional_redundant_cluster_ids={key:[] for key in main_to_redundant_clusters}
main_to_additional_redundant_cluster_descriptions={key:[] for key in main_to_redundant_clusters}

# additional_to_main={key:"" for key in main_to_redundant_clusters}
# print len(cluster_list)
for i in range(1,len(cluster_list)):
    line=cluster_list[i]
    info= line.strip("\n").split("\t")
    curr_accession=info[0]
    curr_cluster_num=info[4][9:-1].strip(" ")
    curr_cluster_id=curr_accession+" "+curr_cluster_num
    curr_redundants1=info[11].strip("\n").split("; ")
    curr_redundants2=info[12].strip("\n").split("; ")
    curr_redundants=curr_redundants1+curr_redundants2
    # print "\nCurrent cluster:",curr_cluster_id

    # If current cluster is considered redundnat to some other main cluster, we record it (and all of its redundant clusters) as redundant to the main cluster,
    # While making sure that these redundant clusters aren't already recorded.
    if curr_cluster_id in redundant_to_main_clusters:

    	# this is the cluster it's redundant to
    	main_clusters=redundant_to_main_clusters[curr_cluster_id]

    	# current cluster can be redundant to several main clusters
    	for each_main_cluster in main_clusters:

	    	# now avoid listing redundant clusters twice

	    	# first, retrieve a list of cluster IDs for clusters that are redundant to main
	    	main_redundant_id_list1=previous_cluster_redundants1[each_main_cluster]
	    	main_redundant_id_list2=previous_cluster_redundants2[each_main_cluster]

	    	# second, retrieve a list of cluster IDS that have already been added as redundant to this main
	    	already_added_redundant_id_list=main_to_additional_redundant_cluster_ids[each_main_cluster]

	    	# third, checking if current cluster hasn't already been added to the list of redundant cluster of this main.
	    	# if it hasn't been then add it to additional redundant list.
	    	if (curr_cluster_id!=each_main_cluster) and (curr_cluster_id not in main_redundant_id_list1) and (curr_cluster_id not in main_redundant_id_list2) and (curr_cluster_id not in already_added_redundant_id_list):
	    		main_to_additional_redundant_cluster_ids[each_main_cluster].append(curr_cluster_id)
	    		main_to_additional_redundant_cluster_descriptions[each_main_cluster].append(short_cluster_descr[curr_cluster_id])
	    		already_added_redundant_id_list=main_to_additional_redundant_cluster_ids[each_main_cluster]


		    # fourth, make a list of redundants of redundant, and for each of them check if they haven't been recorded as redundant to the main cluster yet
	    	for redund_of_red in curr_redundants:
	    		redund_of_red_info=redund_of_red.split(" ")
	    		if len(redund_of_red_info)>2:
		    		redund_of_red_accession=redund_of_red_info[0]
		    		if redund_of_red_info[2]=="Cluster" or redund_of_red_info[2]=="cluster":
						redund_of_red_clusternum=redund_of_red_info[3]
		    		else:
						redund_of_red_clusternum=redund_of_red_info[2][7:]
		    		redund_of_red_cluster_id=redund_of_red_accession+" "+redund_of_red_clusternum

		    		# only record it they haven't been recorded yet
		    		if (redund_of_red_cluster_id!=each_main_cluster) and (redund_of_red_cluster_id not in main_redundant_id_list1) and (redund_of_red_cluster_id not in main_redundant_id_list2) and (redund_of_red_cluster_id not in already_added_redundant_id_list) :
			    		main_to_additional_redundant_cluster_ids[each_main_cluster].append(redund_of_red_cluster_id)
			    		main_to_additional_redundant_cluster_descriptions[each_main_cluster].append(redund_of_red)
			    		already_added_redundant_id_list=main_to_additional_redundant_cluster_ids[each_main_cluster]
	    		else:
					print "short entry in the second loop:",redund_of_red_info


# Changing the table header
header=cluster_list[0].strip("\n")
header+="\tClusters with > 90\% similarity score (chained through transitivity)\n"

# Writing the new file
n_redund_sets=0
n_nonredund=0
with open(output_folder+"/clusterTableNonRedundantNonSimilarNonSequenceSimilar.alldbs.prot_len.correct_domains.txt",'w') as output_file:
    output_file.write(header)
    for i in range(1,len(cluster_list)):
	    line=cluster_list[i]
	    info= line.split("\t")
	    accession=info[0]
	    cluster_num=info[4][9:-1].strip(" ")
	    link=info[5]
	    new_link="http://web.stanford.edu/group/orphan_pks/"+accession+"/index.html#cluster-"+cluster_num
	    info[5]=new_link
	    cluster_id=accession+" "+cluster_num

	    if cluster_id in main_to_additional_redundant_cluster_descriptions:
		    n_redund_sets+=1
		    cluster_description="\t".join(info[0:11])
		    cluster_redundants1=previous_cluster_redundants1[cluster_id]
		    cluster_redundants1_descriptions="; ".join([redundant_cluster_short_descriptions[item_id] for item_id in cluster_redundants1])
		    cluster_redundants2=previous_cluster_redundants2[cluster_id]
		    cluster_redundants2_descriptions="; ".join([redundant_cluster_short_descriptions[item_id] for item_id in cluster_redundants2])
		    output_file.write("\t".join([cluster_description,cluster_redundants1_descriptions,cluster_redundants2_descriptions])+"\t")
		    clusters_to_add=main_to_additional_redundant_cluster_descriptions[cluster_id]
		    for additional_cluster in clusters_to_add:
				output_file.write(additional_cluster+"; ")
		    output_file.write("\n")
	    elif cluster_id in non_redundant_clusters:
		    n_nonredund+=1
		    cluster_description="\t".join(info[0:11])
		    cluster_redundants1=previous_cluster_redundants1[cluster_id]
		    cluster_redundants1_descriptions="; ".join([redundant_cluster_short_descriptions[item_id] for item_id in cluster_redundants1])
		    cluster_redundants2=previous_cluster_redundants2[cluster_id]
		    cluster_redundants2_descriptions="; ".join([redundant_cluster_short_descriptions[item_id] for item_id in cluster_redundants2])
		    output_file.write("\t".join([cluster_description,cluster_redundants1_descriptions,cluster_redundants2_descriptions])+"\t\n")

print "Non-redundants:",n_nonredund
print "Redundant sets:",n_redund_sets



