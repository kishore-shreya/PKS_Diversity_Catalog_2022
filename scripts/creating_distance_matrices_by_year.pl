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


# This script was only slightly modified compared to the initial one.
# It creates the distance matrix for all cluster pairs, including the redundant ones.

# This script is modified so that it uses parsed scores by year as an input instead of 
# pased scores alldbs to avoid unneecssary computation (assuming that parsed scores by year is correct)

use strict;

die "Usage: $0 simFile clusterFile outputfile\n" unless @ARGV >= 3;

my $file = shift @ARGV; # Parsed scores by year
my $clusterFile = shift @ARGV;
my $outputFile = shift @ARGV;

# my @remove_list = ('\#');

my $make_symmetric = 1; ###################################### 

print "Starting to build cluster name array\n";
my @gb_clusternum = build_gb_clusternum($clusterFile); # list of accessions during that year
my $size1=@gb_clusternum;
print "Size of array: $size1\n";

# Build similarities. Accession names as keys otherwise
my %m1_m2_sim = build_m1_m2_sim($file);
my @modules = sort keys %m1_m2_sim;

open(OUT,">$outputFile");

# Print header line
print OUT "Cluster\t";
print OUT join ("\t", @modules);
print OUT "\n";

# Iterate module pairs, print matrix
foreach my $m1 (@modules) {
  print OUT "$m1";
  foreach my $m2 (@modules) {
    my $sim = 0;
    my $sim2 = 0;
    if (exists $m1_m2_sim{$m1}{$m2}){
    	$sim = sprintf("%.3f", $m1_m2_sim{$m1}{$m2});
    } else {
    	print "There's no similarity score for $m1 and $m2 \n";
    }
    if (exists $m1_m2_sim{$m2}{$m1}){
        $sim2 = sprintf("%.3f", $m1_m2_sim{$m2}{$m1});
    } else {
        print "There's no similarity score for $m2 and $m1 \n";
    }
    if($make_symmetric) {
      if($sim2 > $sim) {$sim = $sim2;}
    }
    my $dist = 1-$sim; ###################### Distance, not similarity, required for R hclust function
    print OUT "\t$dist";
  }
  print OUT "\n";
}
close(OUT);
exit;

# Build hash of cluster1 - cluster2 - similarity score
sub build_m1_m2_sim{
  my $file = shift;
  my %m1_m2_sim = ();
  open(IN, $file) || die "Error, can't open file $file\n";
  while(<IN>) {
    chomp;
    my ($m1, $m2, $sim, @rest) = split("\t", $_);
    $m1_m2_sim{$m1}{$m2} = $sim;
  }
  return(%m1_m2_sim);
}
    

sub build_gb_clusternum {
  my $file = shift;
  my @gb_clusternum=();
  open(DESC, $file) || die "Error opening file $file\n";;
  while(<DESC>) {
    chomp;
    my($gb, $clusternum, $date, @rest) = split(" ", $_);
    push (@gb_clusternum,"$gb $clusternum");
  }
  close(DESC);
  return(@gb_clusternum);
}
