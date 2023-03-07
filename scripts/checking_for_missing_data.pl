#!/usr/bin/perl -w

# Copyright (C) 2019 Maureen Hillenmeyer, Aleksandra Nivina

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

# This script is modified from paris to matrix script to find out which pairs are missing from 
# and which pairs are duplicated from the complete parsed score file.

use strict;

die "Usage: $0 simFile clusterFile\n" unless @ARGV >= 3;

my $scoreFile = shift @ARGV;
my $duplicate_outfile = shift @ARGV;
my $missing_outfile = shift @ARGV;

open(DUPE, ">$duplicate_outfile");
open(MISS, ">$missing_outfile");

my %m1_m2_sim = build_m1_m2_sim($scoreFile);
my @modules = sort keys %m1_m2_sim;

my $counter = 0;
foreach my $m1 (@modules) {
	$counter ++;
	foreach my $m2 (@modules) {
        	unless (exists $m1_m2_sim{$m1}{$m2}){
			print MISS "$m1 $m2 \n";
		}
	}
}

print "The num of accessions in score file is: $counter";

close (DUPE);
close (MISS);

exit;


# Build hash of cluster1 - cluster2 - appearance
sub build_m1_m2_sim{
  my $scoreFile = shift;
  my %m1_m2_dict = ();
  open(IN, $scoreFile) || die "Error, can't open file $scoreFile\n";
  while(<IN>) {
    chomp;
    my ($m1, $m2, @rest) = split("\t", $_);
    if (exists $m1_m2_dict{$m1}{$m2}){
	print DUPE "$m1 $m2 \n";
    } else {
	$m1_m2_dict{$m1}{$m2} = 1;
    }
  }
  return(%m1_m2_dict);
}
