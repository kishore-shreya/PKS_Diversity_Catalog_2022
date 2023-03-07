#!/usr/bin/perl -w
use strict;

die "Usage: $0 clusterFile minClusterNum\n" unless @ARGV > 1;

my $clusterFile = shift @ARGV;
my $minClusterNum = shift @ARGV;

open(IN, $clusterFile) || die "Error, can't open file $clusterFile, exiting...\n";
while(<IN>) {
  chomp;
  my($genbank, $desc, $length, $clusterId, $numClusters, $clusterNames, $start, $end) = split("\t", $_);
  next unless $numClusters >= $minClusterNum;
  print "$_\n";
}
