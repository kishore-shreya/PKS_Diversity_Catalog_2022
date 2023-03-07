#!/usr/local/bin/python


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

# This script  that:
# - for the fasta files
# - retrieves the length of amino acid sequence
# - saves these lengths into a new copy of clusterfile

# It takes as arguments:
# - the clusterTable file
# - the location of fasta files
# - the location of output clusterfile

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse


parser = argparse.ArgumentParser(description='Get the range of cluster IDs to parse')

parser.add_argument('-cluster_file','--clusterfile')
parser.add_argument('-inputf','--input_folder')
parser.add_argument('-output_cluster_file','--new_clusterfile')

args=parser.parse_args()

cluster_file=args.clusterfile
my_in_path=args.input_folder
new_cluster_file=args.new_clusterfile


# recording cluster lengths from fasta files, into a dictionary
"""
record_lengths={}
for filename in os.listdir(my_in_path):
    if filename.endswith(".fasta"):
        # print "opening file",my_in_path+filename
        record = SeqIO.read(my_in_path+filename, "fasta")
        record_lengths[record.id]=len(record.seq)
"""

def get_region(cluster_n):
    if len(cluster_n) == 1:
        return 'region00' + cluster_n
    elif len(cluster_n) == 2:
        return 'region0' + cluster_n
    else:
        return 'region' + cluster_n

# loop through lines in the clusterfile
with open(cluster_file,'r') as input_file:
    with open(new_cluster_file,'w') as output_file:
        header=input_file.readline()
        output_file.write(header)
        counter=True
        while counter:
            cluster_record=input_file.readline()
            if cluster_record=="\n" or cluster_record=="":
                counter=False
                break
            else:
                info=cluster_record.split("\t")
                accession=info[0]
                cluster_n=info[4][12:-1]
                region_num = get_region(cluster_n)
                filename = my_in_path + "/" + accession + "." + region_num + ".fasta"
		print "filename", filename
                record = SeqIO.read(filename, "fasta")
                aa_length=len(record.seq)
                # replacing bp_lenth by aa_length
                info[7]=str(aa_length)
                new_record="\t".join(info[:])
                output_file.write(new_record)
