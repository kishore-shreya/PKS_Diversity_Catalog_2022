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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set_style(style="whitegrid")
import matplotlib.colors as mcolors
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
import pandas as pd


known_clusters_main={}
with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/known_clusters/known_clusters/known_clusters_MIBig_and_clusterTable_matched_to_main.txt",'r') as known_cluster_file:
	known_cluster_list=known_cluster_file.readlines()
for line in known_cluster_list:
    cluster_info=line.split("\t")
    cluster_id=cluster_info[0].strip(":")
    product=cluster_info[1].strip("\n")
    known_clusters_main[cluster_id]=0


with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilarNonSequenceSimilar.alldbs.prot_len.correct_domains.txt",'r') as distinct_cluster_file:
	distinct_cluster_list=distinct_cluster_file.readlines()

distinct_cluster_ids={}
for i in range(1,len(distinct_cluster_list)):
	line=distinct_cluster_list[i].split("\t")
	accession=line[0]
	cluster_num=line[4].split(" ")[1][1:-1]
	cluster_id=accession+" "+cluster_num
	distinct_cluster_ids[cluster_id]=0


print "There are",len(known_clusters_main),"known clusters."

known_cluster_indices=[]
distances_to_known_clusters=[]
with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/dataMatrix_alldbs_blastp.txt",'r') as data_matrix_file:
	header=data_matrix_file.readline()
	header_data=header.split("\t")
	# making a list of indices that correspond to known clusters
	for i in range(1,len(header_data)):
		curr_cluster=header_data[i].split(" ")
		curr_accession=curr_cluster[0]
		curr_clusternum=curr_cluster[1]
		curr_cluster_id=curr_accession+" "+curr_clusternum
		if (curr_cluster_id in known_clusters_main) and (curr_cluster_id in distinct_cluster_ids):
			known_cluster_indices.append(i)
			known_clusters_main[curr_cluster_id]=1

	print "There are",len(known_cluster_indices),"known clusters found in the data matrix."
	for key in known_clusters_main:
		if known_clusters_main[key]==0:
			print "Didn't find",key
	end_of_file=False
	while not end_of_file:
		line=data_matrix_file.readline()

		if line=="\n" or line=="":
			end_of_file=False
			break

		line_data=line.split("\t")
		description=line_data[0].split(" ")
		curr_accession=description[0]
		curr_clusternum=description[1]
		curr_cluster_id=curr_accession+" "+curr_clusternum
		# only recording distance info for orphan clusters
		if (curr_cluster_id not in known_clusters_main) and (curr_cluster_id in distinct_cluster_ids):
			max_similarity=0.0
			closest_known=""
			for k in known_cluster_indices:
				similarity=1.0-float(line_data[k])
				if similarity>max_similarity:
					max_similarity=similarity
					closest_known=header_data[k]
			distances_to_known_clusters.append([curr_cluster_id, max_similarity,closest_known])

	print "There are",len(distances_to_known_clusters),"orphan clusters"

distance_df = pd.DataFrame(distances_to_known_clusters, columns=['cluster_id', 'similarity','closest_known'])
distance_df.to_csv("max_similarity_to_known_clusters_distinct.csv")
middle_df = distance_df.loc[(distance_df['similarity']<=0.55)]
middle_df = middle_df.loc[(middle_df['similarity']>=0.45)]
middle_df.to_csv("max_similarity_to_known_distinct_0.45_to_0.55.csv")

fig1 = plt.figure(figsize=(7,6))
ax1 = plt.axes()
pos1 = ax1.get_position()
pos2 = [pos1.x0, pos1.y0 + 0.1,  pos1.width / 2, pos1.height / 1.5]

bins=np.arange(0.0,1.1,0.1)
ax1.set_position(pos2)
ax1.set_xlim(0,1)
ax1.set_xlabel("Max similarity to a known cluster (%)")
ax1.set_xticks(bins)
ax1.set_xticklabels((int(x*100) for x in bins), rotation=90 ) #, rotation=90
# ax1.set_yscale('log')
ax1.set_ylim(0,2500)
ax1.set_ylabel("Number of clusters")
ax1.grid(False)
# ax1.spines['left'].set_color('red')

n, bins, patches=ax1.hist(distance_df['similarity'],bins=bins,color="#5C99D6",align='mid')
print n
print bins
print patches
plt.bar(-10,0,color="#2966A3")
patches[5].set_fc('#2966A3')
patches[6].set_fc('#2966A3')
patches[7].set_fc('#2966A3')
patches[8].set_fc('#2966A3')
patches[9].set_fc('#143352')

percent_truly_orphan=int(np.around(float(n[0]+n[1]+n[2]+n[3]+n[4])/float(len(distances_to_known_clusters))*100,0))
print "percent_truly_orphan",percent_truly_orphan

percent_possibly_homologous=int(np.around(float(n[5]+n[6]+n[7]+n[8])/float(len(distances_to_known_clusters))*100,0))
print "percent_possibly_homologous",percent_possibly_homologous

percent_redundant=int(np.around(float(n[9])/float(len(distances_to_known_clusters))*100,0))
print "percent_redundant",percent_redundant

plt.legend(["Probably diverse ("+str(percent_truly_orphan)+"%)","Possibly homologous ("+str(percent_possibly_homologous)+"%)"],bbox_to_anchor=(pos1.x0+1.1,pos1.y0+0.5,1,0.2))
plt.title("Similarity of distinct orphan clusters to known ones")
plt.savefig("fig_similarity_to_known_clusters_distinct.pdf")

fig2 = plt.figure(figsize=(7,6))
ax1 = plt.axes()
pos1 = ax1.get_position()
pos2 = [pos1.x0, pos1.y0 + 0.1,  pos1.width / 2, pos1.height / 1.5]

ax1.set_position(pos2)
ax1.set_xlim(0,1)
ax1.set_xlabel("Max similarity to a known cluster (%)")
ax1.set_xticks(bins)
ax1.set_xticklabels((int(x*100) for x in bins), rotation=90 ) #, rotation=90
ax1.set_yscale('log')
ax1.set_ylim(1E0,1E4)
ax1.set_ylabel("Number of clusters")
ax1.grid(False)
# ax1.spines['left'].set_color('red')

n, bins, patches=ax1.hist(distance_df['similarity'],bins=bins,color="#5C99D6",align='mid')
plt.bar(-10,1,color="#2966A3")
patches[5].set_fc('#2966A3')
patches[6].set_fc('#2966A3')
patches[7].set_fc('#2966A3')
patches[8].set_fc('#2966A3')
patches[9].set_fc('#143352')
plt.legend(["Probably diverse ("+str(percent_truly_orphan)+"%)","Possibly homologous ("+str(percent_possibly_homologous)+"%)"],bbox_to_anchor=(pos1.x0+0.90,pos1.y0+0.5,1,0.2))
plt.title("Similarity of distinct orphan clusters to known ones")
plt.savefig("fig_similarity_to_known_clusters_log_distinct.pdf")


# percent_truly_orphan 52
# percent_possibly_homologous 48
# percent_redundant 0
