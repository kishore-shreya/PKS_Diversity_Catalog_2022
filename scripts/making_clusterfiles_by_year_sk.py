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

clustertable1994=open("clusterTable1994.nonidentical.txt",'w')
clustertable1995=open("clusterTable1995.nonidentical.txt",'w')
clustertable1996=open("clusterTable1996.nonidentical.txt",'w')
clustertable1997=open("clusterTable1997.nonidentical.txt",'w')
clustertable1998=open("clusterTable1998.nonidentical.txt",'w')
clustertable1999=open("clusterTable1999.nonidentical.txt",'w')
clustertable2000=open("clusterTable2000.nonidentical.txt",'w')
clustertable2001=open("clusterTable2001.nonidentical.txt",'w')
clustertable2002=open("clusterTable2002.nonidentical.txt",'w')
clustertable2003=open("clusterTable2003.nonidentical.txt",'w')
clustertable2004=open("clusterTable2004.nonidentical.txt",'w')
clustertable2005=open("clusterTable2005.nonidentical.txt",'w')
clustertable2006=open("clusterTable2006.nonidentical.txt",'w')
clustertable2007=open("clusterTable2007.nonidentical.txt",'w')
clustertable2008=open("clusterTable2008.nonidentical.txt",'w')
clustertable2009=open("clusterTable2009.nonidentical.txt",'w')
clustertable2010=open("clusterTable2010.nonidentical.txt",'w')
clustertable2011=open("clusterTable2011.nonidentical.txt",'w')
clustertable2012=open("clusterTable2012.nonidentical.txt",'w')
clustertable2013=open("clusterTable2013.nonidentical.txt",'w')
clustertable2014=open("clusterTable2014.nonidentical.txt",'w')
clustertable2015=open("clusterTable2015.nonidentical.txt",'w')
clustertable2016=open("clusterTable2016.nonidentical.txt",'w')
clustertable2017=open("clusterTable2017.nonidentical.txt",'w')
clustertable2018=open("clusterTable2018.nonidentical.txt",'w')
clustertable2019=open("clusterTable2019.nonidentical.txt",'w')
clustertable2020=open("clusterTable2020.nonidentical.txt",'w')
clustertable2021=open("clusterTable2021.nonidentical.txt",'w')
clustertable2022=open("clusterTable2022.nonidentical.txt",'w')

years={1994:clustertable1994,1995:clustertable1995,1996:clustertable1996,1997:clustertable1997,1998:clustertable1998,1999:clustertable1999,\
2000:clustertable2000,2001:clustertable2001,2002:clustertable2002,2003:clustertable2003,2004:clustertable2004,2005:clustertable2005,\
2006:clustertable2006,2007:clustertable2007,2008:clustertable2008,2009:clustertable2009,2010:clustertable2010,2011:clustertable2011,\
2012:clustertable2012,2013:clustertable2013,2014:clustertable2014,2015:clustertable2015,2016:clustertable2016,2017:clustertable2017,\
2018:clustertable2018,2019:clustertable2019,2020:clustertable2020,2021:clustertable2021,2022:clustertable2022}

with open("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/clusterTableNonRedundantNonSimilar.alldbs.prot_len.txt",'r') as clustertable_final:
	counter=True
	line=clustertable_final.readline()
	while counter:
		line=clustertable_final.readline()
		if len(line)>0:		
			dates=[]
			data=line.split("\t")
			try:
				main_date=int(data[2])
			except:
				print "No date for",data[0],"\n"
			else:

				# re-setting the indicators for each clusterFile line
				main_years=[0 for i in range(1994,2023)]
				redundant1_years=[0 for i in range(1994,2023)]
				redundant2_years=[0 for i in range(1994,2023)]
				curr_names=[[] for i in range(1994,2023)]

				# saving the main entry name and date, if it's earlier or equal to a year
				main_entryname=data[0]
				main_date=int(data[2])
				main_clusternum=str(data[4][9:-1])
				for year in years:
					if main_date<=year:
						years[year].write(main_entryname+" "+str(main_clusternum)+" "+str(main_date)+" : "+main_entryname+" "+str(main_clusternum)+" "+str(main_date)+"; ")
						main_years[year-1994]+=1
						curr_names[year-1994].append((main_entryname,main_clusternum))

				# saving the name and date for other entries, if it's earlier or equal to a year
				redundants1=data[19].split("; ")
				redundants2=data[20].split("; ")
				# for entry in redundants1:
				# 	entryname_end=entry.find(" ",1)
				# 	date_start=entry.find("Date=")
				# 	if date_start!=-1:
				# 		entryname=entry[:entryname_end]
				# 		# clusternum_start1=entry.find("cluster",date_start+10,entryname_end+12)
				# 		# clusternum_start2=entry.find("Cluster",date_start+10,entryname_end+12)
				# 		# clusternum_start=np.max(clusternum_start1,clusternum_start2)
				# 		# clusternum_end=entry.find(" ",clusternum_start)
				# 		# clusternum=entry[clusternum_start+7:clusternum_end].strip(" ")

				# 		date=int(entry[date_start+5:date_start+9])
				# 		clusternum=int(entry[date_start+17:entry.find(" ",date_start+18)])

				# 		# loop over years
				# 		for year in years:
				# 			already_present=False
				# 			if date<=year:
				# 				# check if a cluster with a similar name has already been recorded
				# 				for item in curr_names[year-1994]:
				# 					if (entryname.find(item[0])!=-1 and clusternum==item[1]) or (item[0].find(entryname)!=-1 and clusternum==item[1]):
				# 						already_present=True
				# 				# record, if no cluster with a similar name has been recorded
				# 				if not already_present:
				# 					# recording main entry if it hasn't been recorded for this line (to serve as reference later, when we search for distances)
				# 					if main_years[year-1994]==0:
				# 						years[year].write(main_entryname+" "+str(main_clusternum)+" "+str(main_date)+": ")
				# 					years[year].write(entryname+" "+str(clusternum)+" "+str(date)+"; ")
				# 					redundant1_years[year-1994]+=1
				# 					curr_names[year-1994].append((entryname,clusternum))

				for entry in redundants2:
					entryname_end=entry.find(" ",1)
					date_start=entry.find("Date=")
					if date_start!=-1:
						entryname=entry[:entryname_end]
						# clusternum_start1=entry.find("cluster",date_start+10,entryname_end+12)
						# clusternum_start2=entry.find("Cluster",date_start+10,entryname_end+12)
						# clusternum_start=np.max(clusternum_start1,clusternum_start2)
						# clusternum_end=entry.find(" ",clusternum_start)
						# clusternum=entry[clusternum_start+7:clusternum_end].strip(" ")
						date=int(entry[date_start+5:date_start+9])
						clusternum=str(entry[date_start+17:entry.find(" ",date_start+18)])
						clusternum=clusternum[2:-1]
						print "clusternum is", clusternum

						# loop over years
						for year in years:
							already_present=False
							if date<=year:
								# check if a cluster with a similar name has already been recorded
								for item in curr_names[year-1994]:
									if (entryname.find(item[0])!=-1 and clusternum==item[1]) or (item[0].find(entryname)!=-1 and clusternum==item[1]):
										already_present=True
								# record, if no cluster with a similar name has been recorded
								if not already_present:
									# recording main entry if it hasn't been recorded for this line (to serve as reference later, when we search for distances)
									if main_years[year-1994]==0:
										years[year].write(main_entryname+" "+str(main_clusternum)+" "+str(main_date)+" : ")
									years[year].write(entryname+" "+str(clusternum)+" "+str(date)+"; ")
									redundant1_years[year-1994]+=1
									curr_names[year-1994].append((entryname,clusternum))

				# add end of line to the corresponding file if something was recorded for this line of the initial cluaterfile
				for year in years:
					if main_years[year-1994]!=0 or redundant1_years[year-1994]!=0 or redundant2_years[year-1994]!=0:
						years[year].write("\n")
		else:
			break

