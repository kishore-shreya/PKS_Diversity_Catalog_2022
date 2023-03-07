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
import pandas as pd

data=pd.read_csv("/scratch/groups/khosla/Orphan_PKS/results/antismash_results/interesting_clusters/alldbs/Timeline/Similarity/orphan_distinct_cluster_similarities_to_known_timeline_numbers.csv")

labels=[1995,2000,2005,2010,2015,2020]
fig1 = plt.figure(figsize=(7,6))
ax1 = plt.axes()
pos1 = ax1.get_position()
pos2 = [pos1.x0, pos1.y0 + 0.1,  pos1.width / 2, pos1.height / 1.5]
ax1.set_position(pos2)
ax1.set_xlim(1994,2022)
# ax1.set_xlim(1994,2018)
# ax1.set_xlabel("")
ax1.set_xticks(labels)
ax1.set_xticklabels(labels, rotation=90) #, rotation=90
# ax1.set_yscale('log')
# ax1.set_ylim(0,3200)
ax1.set_ylim(0,5000)
ax1.set_ylabel("Number of clusters")
ax1.grid(False)
# ax1.spines['left'].set_color('red')

# plt.plot(years,truly_orphan,color="lightgrey")
# plt.plot(years,possibly_homologs,color="grey")
# plt.plot(years,redundant,color="black")
# plt.plot(years,known,color="red")

plt.plot(data['Year'],data['Truly orphan'],color="lightgrey")
plt.plot(data['Year'],data['Possibly homologous'],color="grey")
plt.plot(data['Year'],data['Redundant'],color="black")
plt.plot(data['Year'],data['Known'],color="red")

plt.legend(["Truly orphan clusters","Possibly homologous\nto known clusters","Redundant\nto known clusters","Known clusters"],bbox_to_anchor=(pos1.x0+1,pos1.y0+0.5,1,0.2))
plt.title("Increase in the number of orphan clusters")
plt.savefig("fig_numbers_of_orphan_clusters_timeline.pdf")

data['Total'] = data['Truly orphan']+data['Possibly homologous']+data['Known']
data['Known_pct'] = data['Known'] / data['Total'] * 100
data['Homologous_pct'] = data['Possibly homologous'] / data['Total'] * 100

print data

labels=[1995,2000,2005,2010,2015,2020]
plt.figure()
fig,ax=plt.subplots(figsize=(6,6))
ax = data.plot(x='Year', y=['Known_pct', 'Homologous_pct'], kind='line', title='Cluster Diversity', legend=False, ylim=(0, 100), color=['red', 'blue'])
plt.xlabel("Year")
plt.ylabel("Clusters (%)")
plt.savefig('fig_cluster_diversity.png')


