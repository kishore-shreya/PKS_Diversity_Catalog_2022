#!/usr/bin/perl -w

use strict;

die "Usage: $0 infolder clusterfile outfile cluster_start_n cluster_end_n\n" unless @ARGV >= 3;

my $redundant = shift;
my $non_redundant = shift;
my $clusterFile = shift;

my %hash = ();
open(LEN, $clusterFile);
while(<LEN>) {
   chomp;
   my(@info) = split("\t", $_);
   my $gb = $info[0];
   my $clusterstring = $info[4]; #print "No value for clusterstring in line $_, file $file";
   $clusterstring =~ /Cluster\s\"(.+)\"/;
   my $clusternum = $1 || $clusterstring;
   my $key_name = "$gb $clusternum";
   $hash{$key_name} = 0;
}
close (LEN);

open(IN, $non_redundant);
while(<IN>) {
  chomp;
  my ($gb1, @rest) = split("\t", $_);
  delete $hash{$gb1};
}
close (IN);

open(IN2, $redundant);
while(<IN2>) {
  chomp;
  my (@accessions) = split("\t", $_);
  foreach my $acc (@accessions){
    delete $hash{$acc};
  }
}
close (IN2);

my $outfile1 = "check_accession_2.txt";
open(OUT1, ">$outfile1");

print OUT1 "Printing accessions that didn't appear in both files";
foreach my $no_match (keys %hash){
    print OUT1 "No result: $no_match\n";
}
