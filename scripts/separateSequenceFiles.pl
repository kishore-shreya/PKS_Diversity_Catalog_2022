#!/usr/bin/perl -w

# Copyright (C) 2019 Maureen Hillenmeyer

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

use strict;

use Bio::SeqIO;

die "Usage: $0 infile format outfolder\n" unless @ARGV >2;

my $file = shift @ARGV;
my $format = shift @ARGV;
my $outfolder = shift @ARGV;
my $outformat = $format;
if($format eq 'genbank') {$outformat = "gb";}

my $seqin = Bio::SeqIO->new(-file   => $file,
                            -format => $format);
 
`mkdir $outfolder`;
while (my $seq = $seqin->next_seq) {

  my $id = $seq->id;
  #$id =~ /ref\|(.*)\|/;
  #$id =~ /gb\|(.*)\|/;
  #$id =~ /dbj\|(.*)\|/;
  $id =~ /gi\|(.*?)\|/;
  my $accession = $1 || $id;
  warn "id = $id, accession = $accession\n";
  
  my $seqout = Bio::SeqIO->new(-file   => ">$outfolder/$accession.$outformat",
                            -format => $format);
 
  $seqout->write_seq($seq);
}