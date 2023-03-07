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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
import matplotlib.cm as cm
import scipy.cluster.hierarchy
import argparse
# import fastcluster

# Plotting all clusters
# python plot_dendrograms.py -show_all Yes -show_known Yes -only_show_known No -show_accession Yes -show_species Yes -show_date Yes -show_type Yes -show_domains Yes -color_known Yes > dendrogram_output.txt

# Plotting only unique clusters
# python plot_dendrograms.py -show_all No -show_known Yes -only_show_known No -show_accession Yes -show_species Yes -show_date Yes -show_type Yes -show_domains Yes -color_known Yes > dendrogram_output_only_known.txt

parser = argparse.ArgumentParser(description='What to display')

# Parsing the command-line argument
parser.add_argument('-show_all','--showall')
parser.add_argument('-show_known','--showknown')
parser.add_argument('-only_show_known','--onlyshowknown')

parser.add_argument('-show_accession','--showaccession')
parser.add_argument('-show_species','--showspecies')
parser.add_argument('-show_date','--showdate')
parser.add_argument('-show_type','--showtype')
parser.add_argument('-show_domains','--showdomains')
parser.add_argument('-color_known','--colorknown')

args=parser.parse_args()

show_all=args.showall
show_known=args.showknown
only_show_known=args.onlyshowknown

show_accession=args.showaccession
show_species=args.showspecies
show_date=args.showdate
show_type=args.showtype
show_domains=args.showdomains
color_known=args.colorknown


# ##########################################
# ##### Extracting redundant clusters ######
if show_all=="No":
    main_redundant_clusters={}
    other_redundant_clusters={}
    nonredundant_clusters={}

    with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs.redundant_chained_new.cutoff0.9.txt",'r') as redundancy_file:
        redundancy_list=redundancy_file.readlines()

    for k1 in range(0,len(redundancy_list)):

        # Splitting the line into a list, where each element is either the main cluster (first element of the list)
        # or a redundant cluster (other elements in the list)
        redundant_set=redundancy_list[k1].split("\t")
        redundant_main=redundant_set[0]
        main_redundant_clusters[redundant_main]=1
        for i in range(1,len(redundant_set)):
            redundant_other=redundant_set[i]
            other_redundant_clusters[redundant_other]=1

    with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs.nonredundant_chained_new.cutoff0.9.txt",'r') as nonredundancy_file:
        nonredundancy_list=nonredundancy_file.readlines()

    for k1 in range(0,len(nonredundancy_list)):

        # Non-redundant clusters
        nonredundant_cluster=nonredundancy_list[k1].strip("\n")
        nonredundant_clusters[nonredundant_cluster]=1


###########################################################################################################################
##### Going through the correct distance matrix file and extracting accession and cluster number, to match to a description #######
cluster_descriptions={}
with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.correct_domains.txt",'r') as cluster_table_file:
    cluster_table=cluster_table_file.readlines()
    cluster_table=cluster_table[1:]
    for line in cluster_table:
        line_info=line.split("\t")
        accession=line_info[0]
        date=int(line_info[2])
        species=line_info[3]
        clusternum=line_info[4][9:-1]
        n_domains=line_info[9]
        clus_type=line_info[10]
        cluster_descriptions[(accession,clusternum)]=(date,species,clus_type,n_domains)


# ######################################
# ##### Extracting known clusters ######
known_clusters_all={}
known_clusters_main={}

if show_known=="Yes":

    if show_all=="Yes":   

        with open("known_clusters_MIBig_and_clusterTable_matched_to_all_leaves.txt",'r') as known_file:
            known_cluster_list=known_file.readlines()
            for line in known_cluster_list:
                cluster_info=line.split("\t")
                cluster_id=cluster_info[0].strip(":")
		product=cluster_info[1].strip("\n")
                known_clusters_all[cluster_id]=product

    elif show_all=="No":

        with open("known_clusters_MIBig_and_clusterTable_matched_to_main.txt",'r') as known_file:
            known_cluster_list=known_file.readlines()
            for line in known_cluster_list:
                cluster_info=line.split("\t")
                cluster_id=cluster_info[0].strip(":")
                product=cluster_info[1].strip("\n")
                known_clusters_main[cluster_id]=product

###########################################################################################################################
##### Going through the distance matrix file and extracting accession and cluster number, to match to a description #######

with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/dataMatrix_alldbs_blastp.txt",'r') as distance_file:
    distance_data=distance_file.readlines()
tbd=[] # to-be-deleted
distance_matrix=[]
id_labels=[]
labels=[]
counter=0

# First line is column names (clusters).

# Going through the rest of the data matrix records
for k in range(1,len(distance_data)):
    line=distance_data[k]
    line=line.strip('\n').split("\t")
    
    # distances (excluding the first column which is descriptive) are kept in the distance_matrix
    distances=line[1:]
    distance_matrix.append(distances)
    
    # the rest of information is split
    description=line[0].split(" ")
    accession=description[0]
    clusternum=description[1]
    cluster_id=accession+" "+clusternum
    cluster_full_description=cluster_descriptions[(accession,clusternum)]
    date=cluster_full_description[0]
    species=cluster_full_description[1]
    clus_type=cluster_full_description[2]
    n_domains=cluster_full_description[3]
    

    curr_description=""
    if show_accession=="Yes":
        curr_description+=accession+" "+clusternum+" "
    if show_species=="Yes":
        curr_description+=species+" "
    if show_date=="Yes":
        curr_description+=str(date)+" "
    if show_type=="Yes":
        if len(clus_type) >= 100:
            print "Omitting type for accession: ",accession," cluster: ",clusternum
            curr_description+=clus_type[0:100]+" (omitted) "
        else:
            curr_description+=clus_type+" "
    if show_domains=="Yes":
        curr_description+=n_domains
       

    if show_known=="Yes":

        if show_all=="Yes":
            if cluster_id in known_clusters_all:
               curr_description+=", "+known_clusters_all[cluster_id]+" synthase"
        elif show_all=="No":
            if cluster_id in known_clusters_main:
               curr_description+=", "+known_clusters_main[cluster_id]+" synthase"

    if show_all=="No":
        if cluster_id in main_redundant_clusters:
           curr_description+=", seen more than once"
        if (cluster_id not in main_redundant_clusters) and (cluster_id not in nonredundant_clusters):
            tbd.append(counter) 

    if only_show_known=="Yes":
        if cluster_id not in known_clusters_all:
            tbd.append(counter)

    counter+=1
    labels.append(unicode(curr_description, "utf-8"))
    id_labels.append(accession+" "+clusternum)
        

#################################################################################
##### Deleting distances and labels that correspond to redundant clusters #######

# transforming data into arrays
matrix = np.array(distance_matrix,ndmin=2).astype(np.float)
labels=np.array(labels,ndmin=1)
id_labels=np.array(id_labels,ndmin=1)

if show_all=="No" or only_show_known=="Yes":

    # deleting rows that correspond to redundant clusters
    matrix = np.delete(matrix,tbd,axis=0)

    # deleting labels that correspond to redundant clusters
    labels=np.delete(labels,tbd,axis=0)
    id_labels=np.delete(id_labels,tbd,axis=0)
    print "only",len(labels),"leaves left"

    # deleting columns that correspond to redundant clusters
    matrix = np.delete(matrix,tbd,axis=1)

# diagonal is nul
np.fill_diagonal(matrix,0)

#####################################
##### Plotting the dendrogram #######

if only_show_known=="Yes":
    interesting_labels=["KU568466 r1c1","AY118081 r1c1","AB089954 r1c1","AY509120 r1c1","MF033535 r1c1","KP997155 r1c1","KP997155 r1c1","AR897801 r1c1","EU220288 r1c1","AF016585 r1c1"]
    interesting_labels_full=["KU568466 cluster 1","AY118081 cluster 1","AB089954 cluster 1","AY509120 cluster 1","MF033535 cluster 1","KP997155 cluster 1","KP997155 cluster 1","AR897801 cluster 1","EU220288 cluster 1","AF016585 cluster 1"]
    interesting_indices=[]
    for k in range(0,len(id_labels)):
        if id_labels[k] in interesting_labels:
            interesting_indices.append(k)

elif color_known=="Yes":
    known_labels_full=[]
    for item in id_labels:
        if show_all=="Yes":
            if item in known_clusters_all:
                item_data=item.split(" ")
                full_label=item_data[0]+" "+item_data[1]
                known_labels_full.append(full_label)
        else:
            if item in known_clusters_main:
                item_data=item.split(" ")
                full_label=item_data[0]+" "+item_data[1]
                known_labels_full.append(full_label)


# df = pd.DataFrame(data=matrix,columns = ["attr_%d" % j for j in range(matrix.shape[1])])

matplotlib.rcParams['lines.linewidth'] = 0.3

Z = linkage(scipy.spatial.distance.squareform(matrix),method='weighted')

if only_show_known=="Yes":
    fig,axes = pyplot.subplots(figsize=(10,12))
elif show_all=="Yes":
    fig,axes = pyplot.subplots(figsize=(15,200))
else:
    fig,axes = pyplot.subplots(figsize=(15,120))

dflt_col = "grey" 
leaf_colors={}

if only_show_known=="Yes":
    for i in range(0,len(id_labels)):
        if i in interesting_indices:
            leaf_colors["attr_"+str(i)]="red"
        else:
            leaf_colors["attr_"+str(i)]=dflt_col
        if i in interesting_indices:
            leaf_colors["attr_"+str(i)]="red"
        else:
            leaf_colors["attr_"+str(i)]=dflt_col

else:
    for i in range(0,len(id_labels)):
        leaf_colors["attr_"+str(i)]=dflt_col



link_cols = {}
for i, i12 in enumerate(Z[:,:2].astype(int)):
    # print i
    # print i12
    # for x in i12:
    #     if x<=len(Z):
    #         print x, leaf_colors["attr_%d"%x]
    #     else:
    #         print x, "larger than Z"
    c1, c2 = (link_cols[x] if x > len(Z) else leaf_colors["attr_%d"%x] for x in i12)
    # print c1, c2
    link_cols[i+1+len(Z)] = c1 if c1 == c2 else dflt_col
    # print "link color for ",i+1+len(Z),link_cols[i+1+len(Z)]
    # print "\n"

dg =dendrogram(Z, truncate_mode=None, orientation='right',
    #leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=1,  # font size for the x axis labels
    labels = labels,
    color_threshold=0,
    link_color_func=lambda x: link_cols[x],
    ax=axes
)

if only_show_known=="Yes":
    axes.tick_params(axis='y', which='major', labelsize=2)        
elif show_all=="Yes":
    axes.tick_params(axis='y', which='major', labelsize=0.5)
else:
    axes.tick_params(axis='y', which='major', labelsize=1)

y_labels = axes.get_ymajorticklabels()

if only_show_known=="Yes":
    for label in y_labels:
        label.set_color('black')
        if any(s in label.get_text() for s in interesting_labels_full):
            label.set_color('red')
elif color_known=="Yes":
    for label in y_labels:
        label.set_color('black')
        if any(s in label.get_text() for s in known_labels_full):
            label.set_color('red')    

pyplot.title('Hierarchical Clustering Dendrogram')
axes.set_xlabel('Distance Score blastp')

# If you need a horizontal line, we uncommend the following
#if show_all=="Yes":
#    axes2=axes.twinx()
#    axes2.set_ylim(0,200)
#    axes2.plot([0.1,0.1],[0,190],color="grey",linewidth=0.6)

fig.tight_layout() #rect=[0, 0.03, 0.95, 1]

name="dendrogram"
if show_all=="No":
    name+="_distinct"
else:
    name+="_all"
if only_show_known=="Yes":
    name+="_onlyknown"
elif show_known=="Yes":
    name+="_known"
if show_accession=="Yes":
    name+="_accession"
if show_species=="Yes":
   name+="_species"
if show_date=="Yes":
    name+="_date"
if show_type=="Yes":
    name+="_type"
if show_domains=="Yes":
    name+="_domains"

pyplot.savefig(name+'.pdf') 

