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


parsedscores1994=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.1994",'w')
parsedscores1995=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.1995",'w')
parsedscores1996=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.1996",'w')
parsedscores1997=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.1997",'w')
parsedscores1998=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.1998",'w')
parsedscores1999=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.1999",'w')
parsedscores2000=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2000",'w')
parsedscores2001=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2001",'w')
parsedscores2002=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2002",'w')
parsedscores2003=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2003",'w')
parsedscores2004=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2004",'w')
parsedscores2005=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2005",'w')
parsedscores2006=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2006",'w')
parsedscores2007=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2007",'w')
parsedscores2008=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2008",'w')
parsedscores2009=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2009",'w')
parsedscores2010=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2010",'w')
parsedscores2011=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2011",'w')
parsedscores2012=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2012",'w')
parsedscores2013=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2013",'w')
parsedscores2014=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2014",'w')
parsedscores2015=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2015",'w')
parsedscores2016=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2016",'w')
parsedscores2017=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2017",'w')
parsedscores2018=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2018",'w')
parsedscores2019=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2019",'w')
parsedscores2020=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2020",'w')
parsedscores2021=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2021",'w')
parsedscores2022=open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/ParsedScores/parsed_scores.2022",'w')

parsed_score_files={1994:parsedscores1994,1995:parsedscores1995,1996:parsedscores1996,1997:parsedscores1997,1998:parsedscores1998,1999:parsedscores1999,\
2000:parsedscores2000,2001:parsedscores2001,2002:parsedscores2002,2003:parsedscores2003,2004:parsedscores2004,2005:parsedscores2005,\
2006:parsedscores2006,2007:parsedscores2007,2008:parsedscores2008,2009:parsedscores2009,2010:parsedscores2010,2011:parsedscores2011,\
2012:parsedscores2012,2013:parsedscores2013,2014:parsedscores2014,2015:parsedscores2015,2016:parsedscores2016,2017:parsedscores2017,\
2018:parsedscores2018,2019:parsedscores2019,2020:parsedscores2020,2021:parsedscores2021,2022:parsedscores2022}

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
			for p in range(0,len(cluster_years)):
				year=cluster_years[p]
				if cluster1 in clusters[p] and cluster2 in clusters[p]:
					parsed_score_files[year].write("\t".join([cluster1,cluster2,similarity]))
					parsed_score_files[year].write("\n")






