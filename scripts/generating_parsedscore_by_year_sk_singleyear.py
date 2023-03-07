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


clusters=[]
cluster_year=2016


with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ClusterTablesTimeline/clusterTable"+str(cluster_year)+".nonidentical.txt",'r') as cluster_table:
	file_end=False
	while not file_end:
		line=cluster_table.readline()
		if line=="\n" or line=="":
			file_end=True
		else:
			line=line.split(" : ")
			main_cluster=line[0].split(" ")
			main_accession=main_cluster[0]
			main_clusternum=main_cluster[1]
			clusters.append(main_accession+" "+main_clusternum)

scorefile_end=False

PATH = "/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/singleyear_parsed_scores." + str(cluster_year)
parsed_score_files = {cluster_year: open(PATH, "w")}

with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores_alldbs") as score_file:
	while not scorefile_end:
		line=score_file.readline()
		if line=="\n" or line=="":
			scorefile_end=True
		else:
			info=line.split("\t")
			cluster1=info[0]
			cluster2=info[1]
			similarity=str(round(float(info[2]),3))
			if cluster1 in clusters and cluster2 in clusters:
				parsed_score_files[cluster_year].write("\t".join([cluster1,cluster2,similarity]))
				parsed_score_files[cluster_year].write("\n")






