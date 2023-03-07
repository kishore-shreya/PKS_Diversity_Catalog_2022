#!/usr/local/bin/python3.7


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

import json
import os
import csv
import numpy as np
directory = "./mibig_json_1.4/"



cluster_table=[]

for filename in os.listdir(directory):
    with open(os.path.join(directory, filename)) as f:
        # print (os.path.join(directory, filename))
        data = json.load(f)
    
    # There are 2 ways how a PKS is annotated.
    # Sometimes both are annotated as Polyketide, but sometimes only one of them is,
    # so we have to consider all cases
    # c = data.get('general_params').get('Polyketide') # case #1: the product is annotated as "Polyketide" (c!='None')
    MIBig_num=data.get('general_params').get("mibig_accession")
    biosyn_class=data['general_params']['biosyn_class'][0]
    i = data.get('general_params').get('loci').get('nucl_acc')
    accession=i[0].get('Accession')
    compound = data.get('general_params').get('compounds')[0].get('compound')
    cluster_table.append([accession,MIBig_num,compound]) # case #2: the biosynthetic class includes "Polyketide" ('Polyketide' in biosyn_class)
    
 

# saving other PKS in the second file: we might find some of them in the dendrogram, but not necessarily (if they're not assembly-line PKS)
np.savetxt("all_clusters_mibig.txt", cluster_table, delimiter=",", fmt='%s')



