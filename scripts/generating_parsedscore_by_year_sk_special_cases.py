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


clusters=[[] for i in range(29)]
cluster_years=range(1994,2023)
# cluster_years=range(1994,1996)


for k in range(0,len(cluster_years)):
	year=cluster_years[k]
	# print year
	with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ClusterTablesTimeline/clusterTable"+str(year)+".nonidentical.txt",'r') as cluster_table:
		file_end=False
		while not file_end:
			line=cluster_table.readline()
			if line=="\n" or line=="":
				file_end=True
			else:
				line=line.split(" : ")
				# print line
				main_cluster=line[0].split(" ")
				# print main_cluster
				main_accession=main_cluster[0]
				main_clusternum=main_cluster[1]
				# print clusters[k]
				# print "\n"
				clusters[k].append(main_accession+" "+main_clusternum)

print clusters[0]
print clusters[1]

scorefile_end=False

PATH = "/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/special_parsed_scores."
parsed_score_files = {}
for year in cluster_years:
	this_path = PATH + str(year)
	parsed_score_files[year] = open(this_path, "w")

with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/prot_fasta_withClusterN_for_blastp_alldbs/parsed_scores/parsed_scores.special_cases") as score_file:
	while not scorefile_end:
		line=score_file.readline()
		if line=="\n" or line=="":
			scorefile_end=True
		else:
			info=line.split("\t")
			cluster1=info[0]
			cluster2=info[1]
			similarity=str(round(float(info[2]),3))
			for p in range(0,len(cluster_years)):
				year=cluster_years[p]
				if cluster1 in clusters[p] and cluster2 in clusters[p]:
					parsed_score_files[year].write("\t".join([cluster1,cluster2,similarity]))
					parsed_score_files[year].write("\n")






