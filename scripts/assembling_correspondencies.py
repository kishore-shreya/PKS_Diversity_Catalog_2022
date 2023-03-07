
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

with open("known_all_clusters_from_MIBig.txt",'r') as MIBig_cluster_file:
	MIBig_clusters=MIBig_cluster_file.readlines()
with open("known_PKS_clusters_from_clusterTable_edited.txt",'r') as table_cluster_file:
	table_clusters=table_cluster_file.readlines()

known_clusters={}
for line in MIBig_clusters:
	known_info=line.split("\t")
	cluster_id=known_info[0]
	product=known_info[1].strip("\n")
	if cluster_id in known_clusters:
		prev_products=known_clusters[cluster_id].split("; ")
		for prev_product in prev_products:
			if prev_product.find(product)!=-1 or product.find(prev_product)!=-1:
				pass
			else:
				prev_products.append(product)
		new_products="; ".join(prev_products)
		known_clusters[cluster_id]=new_products
	else:
		known_clusters[cluster_id]=product

for line in table_clusters:
	known_info=line.split("\t")
	cluster_id=known_info[0]
	product=known_info[1].strip("\n")
	if cluster_id in known_clusters:
		prev_products=known_clusters[cluster_id].split("; ")
		for prev_product in prev_products:
			if prev_product.find(product)!=-1 or product.find(prev_product)!=-1:
				pass
			else:
				prev_products.append(product)
		new_products="; ".join(prev_products)
		known_clusters[cluster_id]=new_products
	else:
		known_clusters[cluster_id]=product


with open("known_clusters_MIBig_and_clusterTable.txt",'w') as output_file:
	for key in known_clusters:
		output_file.write(key+"\t")
		products=known_clusters[key].split("\t")
		product_desc=products[0]
		if len(products)>1:
			for i in range(0,len(products)):
				product_desc+=" and "+products[i]
		output_file.write(product_desc+"\n")

