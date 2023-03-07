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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

redund_file_list = []
nonredund_file_list = []
for year in range(1994, 2023):
	redund_path = "/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores." + str(year) + ".redundant_chained.cutoff0.9.txt"
	nonredund_path = "/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores." + str(year) + ".nonredundant_chained.cutoff0.9.txt"
	redund_file_list.append(redund_path)
	nonredund_file_list.append(nonredund_path)

# Get all known clusters
lst_known_cluster = []
with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/known_clusters/known_clusters/known_clusters_MIBig_and_clusterTable_edited.txt",'r') as known_cluster_file:
	known_cluster_list=known_cluster_file.readlines()
for line in known_cluster_list:
    cluster_info=line.split("\t")
    cluster_id=cluster_info[0].strip(":")
    lst_known_cluster.append(cluster_id)
known_cluster_file.close()

with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/known_clusters/known_clusters/known_clusters_MIBig_and_clusterTable_matched_to_all_leaves.txt",'r') as known_cluster_file2:
	known_cluster_list2=known_cluster_file2.readlines()
for line in known_cluster_list2:
    cluster_info=line.split("\t")
    cluster_id=cluster_info[0]
    if cluster_id not in lst_known_cluster:
	    lst_known_cluster.append(cluster_id)
known_cluster_file2.close()

with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/known_clusters/known_clusters/known_clusters_MIBig_and_clusterTable_matched_to_main.txt",'r') as known_cluster_file3:
	known_cluster_list3=known_cluster_file3.readlines()
for line in known_cluster_list3:
    cluster_info=line.split("\t")
    cluster_id=cluster_info[0]
    if cluster_id not in lst_known_cluster:
	    lst_known_cluster.append(cluster_id)
known_cluster_file3.close()	  

non_red_numbers=[0 for i in range(1994,2023)]
main_red_numbers=[0 for i in range(1994,2023)]
red_numbers = [0 for i in range(1994, 2023)]
known_uniques = [0 for i in range(1994, 2023)]

print len(lst_known_cluster)

print lst_known_cluster

years=range(1994,2023)
for i in range(0, len(years)):

	known = 0

	nonredund_file=open(nonredund_file_list[i],'r')
	nonredund_clusters=[non_red.strip("\n") for non_red in nonredund_file.readlines() if non_red.strip()]
	nonredund_file.close()
	non_red_numbers[i]=len(nonredund_clusters)

	redund_file=open(redund_file_list[i],'r')
	redund_cluster_sets=redund_file.readlines()
	redund_file.close()

	main_redundant_clusters = []
	lst_redundant_clusters = []
	for redundant_set in redund_cluster_sets:
		redundant_list=redundant_set.strip("\n").split("\t")
		if redundant_list:
			main_redundant_clusters.append(redundant_list[0].strip(":").strip("\t"))
			lst_redundant_clusters += redundant_list[1:]
	main_red_numbers[i] = len(redund_cluster_sets)
	lst_redundant_clusters = list(set(lst_redundant_clusters))
	red_numbers[i] = len(lst_redundant_clusters)

	for known_cluster in lst_known_cluster:
		if known_cluster in nonredund_clusters or known_cluster in main_redundant_clusters:
			known += 1

	known_uniques[i] = known

dic = {"Year": list(years), "NonRed": non_red_numbers, "MainRed": main_red_numbers, "Red": red_numbers, "Known": known_uniques}
data = pd.DataFrame(dic)
data['Unique'] = data['NonRed'] + data['MainRed']
data['Total'] = data['Unique'] + data['Red']
data['Unknown'] = data['Unique'] - data['Known']
data['Rediscovery Rate'] = 100*data['Red'] / data['Total']

# Columns = ['Year', 'NonRed', 'MainRed', 'Red', 'Known', 'Unique', 'Total', 'Unknown']
data.to_csv("timeline_unique_red.csv")

# Plotting the number of unique clusters
plt.figure()
fig,ax=plt.subplots(figsize=(6,6))
ax = data.plot(x='Year', y='Unique', kind='line', title='Unique clusters by year', legend=False, ylim=(0, 10000))
plt.xlabel("Year")
plt.ylabel("The number of unique clusters")
plt.savefig('fig_unique_cluster_by_year.png')

# Plotting the known clusters vs. unique clusters in log scale
plt.figure()
fig,ax=plt.subplots(figsize=(6, 6))
ax = data.plot.bar(x='Year', y=['Known', 'Unknown'], stacked=True, log=True, ylim=(0.5, 1E4), color=['red', 'blue'], title='Known clusters in unique clusters by year (log scale)')
plt.ylabel("The number of known clusters in unique clusters")
plt.xlabel("Year")
plt.savefig("fig_known_by_year.png")

# Plotting the known clusters vs. unique clusters in log scale
plt.figure()
fig,ax=plt.subplots(figsize=(6, 6))
ax = data.plot.bar(x='Year', y=['Known', 'Unknown'], stacked=True, ylim=(0, 1E4), color=['red', 'blue'], title='Known clusters in unique clusters by year')
plt.ylabel("The number of known clusters in unique clusters")
plt.xlabel("Year")
plt.savefig("fig_known_by_year_non_log_scale.png")


# Plotting the percentage of unique clusters among total clusters
plt.figure()
fig,ax=plt.subplots(figsize=(6, 6))
ax = data.plot(x='Year', y='Rediscovery Rate')
plt.xlabel("Year")
plt.ylabel("Redundant clusters / total clusters")
plt.savefig("fig_red_pct_by_year.png")
