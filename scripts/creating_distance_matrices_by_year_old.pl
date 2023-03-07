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

use strict;

die "Usage: $0 simFile clusterFile outputfile\n" unless @ARGV >= 2;

my $file = shift @ARGV; # "moduleSimilarity_0-478.txt";
my $clusterFile = shift @ARGV;
my $outputFile = shift @ARGV;
# my $selected_file = 0;
# if($clusterFile =~ /selected/) {$selected_file = 1;} # Format the labels differently for hand-selected clusters

# my @remove_list = ('\#');

my $make_symmetric = 1; ###################################### 

# my $gene_cluster_phrase = "gene cluster";
# my $gene_cluster_notation = "*** GENE CLUSTER ***";
#my $gene_cluster_notation = "";
# my $desc_length = 120;

print "Starting to build cluster name array\n";
my @gb_clusternum = build_gb_clusternum($clusterFile); # key = gb, key2 = clusternum, val = desc
my $size1=@gb_clusternum;
print "Size of array: $size1\n";
# my %gb_clusternum_desc = %$gb_clusternum_desc;
# my %gb_clusternum_date = %$gb_clusternum_date;
# my %gb_clusternum_domains = %$gb_clusternum_domains;

# old : gb only
#my %gb_known_gene_cluster = build_gb_known_gene_cluster($clusterFile); # Find which gb descs have the term "gene cluster" in it
# new gb and clusternum
# my %gb_clusternum_known_gene_cluster = build_gb_known_gene_cluster($clusterFile); # Find which gb descs have the term "gene cluster" in it



# Build similarities.  This is where the name/description is made
my %m1_m2_sim = build_m1_m2_sim($file,@gb_clusternum);
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
    my $sim =sprintf("%.3f", $m1_m2_sim{$m1}{$m2}) || 0;
    my $sim2 = sprintf("%.3f", $m1_m2_sim{$m2}{$m1}) || 0;
    #if($m1 =~ "AY661566" && $m2 =~ "X62569") { ################################ debug
    #  die "sim = $sim, sim2 = $sim2\n";
    #}
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
  my @gb_clusternum=@_;
  my $length2=@gb_clusternum;
  print "Length of array in subroutine: $length2\n";
  my %m1_m2_sim = ();
  open(IN, $file) || die "Error, can't open file $file\n";
  #my $header = <IN>;
  while(<IN>) {
    chomp;
    my ($m1, $m2, $sim, @rest) = split("\t", $_);

    # my($gb1, $clusternum1) = split(' ', $m1);
    # my($gb2, $clusternum2) = split(' ', $m2);
    # my $label1="gb1 clusternum1";
    # my $label2="gb2 clusternum2";
    my $match1=0;
    my $match2=0;
    for my $item (@gb_clusternum) {
      if ($m1 eq $item) {
        $match1=1;
      }
      if ($m2 eq $item) {
        $match2=1;
      }
    }
    if ($match1 and $match2) {


    # my $desc1 = $gb_clusternum_desc{$gb1}{$clusternum1};
    # my $desc2 = $gb_clusternum_desc{$gb2}{$clusternum2};
    # my $date1 = $gb_clusternum_date{$gb1}{$clusternum1};
    # my $date2 = $gb_clusternum_date{$gb2}{$clusternum2};
    # my $domains1 = $gb_clusternum_domains{$gb1}{$clusternum1};
    # my $domains2 = $gb_clusternum_domains{$gb2}{$clusternum2};
    #warn "$gb1 $clusternum1 $desc1\t$gb2 $clusternum2 $desc2\n";

    # foreach my $remove(@remove_list){
    #   $desc1 =~ s/$remove//g;
    #   $desc2 =~ s/$remove//g;
    # }


    ##############################################
    # Labels
    # 
    # my $label1 = "$m1 $date1 $domains1 $desc1";
    # my $label2 = "$m2 $date2 $domains2 $desc2";

    # Labels for clusters labeled "gene clsuter"
    ############### TODO for known ones - change this to have the known number of domains / date?
    #if(my $known1 = $gb_known_gene_cluster{$gb1}) {
    #if(my $known1 = $gb_known_gene_cluster{$gb1}) {
    # if(my $known1 = $gb_clusternum_known_gene_cluster{$gb1}{$clusternum1}) {
    #   #$label1 = "*** $known1 *** $m1$desc1";
    #   $label1 = "$gene_cluster_notation $m1 $date1 $domains1 $known1";# $desc1";
    #   if($known1 ne $desc1) {$label1 .= " $desc1";}
    # }
    #if(my $known2 = $gb_known_gene_cluster{$gb2}) {
    #if(my $known2 = $gb_known_gene_cluster{$gb2}) {
    # if(my $known2 = $gb_clusternum_known_gene_cluster{$gb2}{$clusternum2}) {
    #   #$label2 = "*** $known2 *** $m2$desc2";
    #   $label2 = "$gene_cluster_notation $m2 $date2 $domains2 $known2";# $desc2";
    #   if($known2 ne $desc2) {$label2 .= " $desc2";}
    # }

    #########################HERE
    # if($selected_file) {
    #   my($name1, $rest1) = split("\t", $desc1);
    #   my($name2, $rest2) = split("\t", $desc2);
    #   $label1 = "$name1 ($m1) $rest1";
    #   $label2 = "$name2 ($m2) $rest2";
    # }


      $m1_m2_sim{$m1}{$m2} = $sim;
     }
  }
  return(%m1_m2_sim);
}
    

sub build_gb_clusternum {
  my $file = shift;
  my @gb_clusternum=();
  # my %clusternum =  ();
  # my %dates =  ();
  # my %domains =  ();
  open(DESC, $file) || die "Error opening file $file\n";;
  while(<DESC>) {
    chomp;
    my($gb, $clusternum, $date, @rest) = split(" ", $_);
    push (@gb_clusternum,"$gb $clusternum");
    # $date =~ tr/;\s//d
    # my $cluster_phrase= $rest[1];
    # $cluster_phrase =~ /Cluster (\d+)/ || warn "Warning, clsuter phrase = $cluster_phrase\n";
    # my $clusternum = $1;
    # my $ks = $rest[6] || "";#warn "Warning, no KS in $file, line = $_\n";
    # my $at = $rest[7] || "";#warn "Warning, no AT in $file, line = $_\n";
    # my $c = $rest[12];
    # my $a = $rest[13];

    # if($selected_file) {
    #   my($name, @rest) = split(' ', $desc);
    #   #warn "name = $name";
    #   #my $rest = join (' ', @rest);
    #   my $genus = shift @rest;
    #   my $species = shift @rest;
    #   my $species2 = "";
    #   if($species eq "sp.") {$species2 = shift @rest;}
    #   if($desc =~ /atent (.*)/) {$genus = "Patent"; $species = $1;}
    #   #$desc = "$name\t$date, $ks KS, $at AT, $c C, $a A, $rest";
    #   $desc = "$name\t$date, $ks KS, $at AT, $c C, $a A, $genus $species $species2";
      
    # }
    # else {
    #   #$desc = ", $date, $ks KS, $at AT, $c C, $a A, $desc)";
      
    # }

    # my $domains = "$ks KS, $at AT, $c C, $a A";

    # if((length $desc) > $desc_length) {
    #   $desc = substr($desc, 0, $desc_length);
    #   $desc = $desc."...";
    # }
    # $desc =~ s/\'//g;

    #warn "gb = $gb, desc = $desc\n";
    # $desc{$gb}{$clusternum} = $desc;
    # $dates{$gb}{$clusternum} = $date;
    # $domains{$gb}{$clusternum} = $domains;
  }
  close(DESC);
  return(@gb_clusternum);
}

# Oct 2013 edits...
# sub build_gb_known_gene_cluster{
#   my $file = shift;
#   my %hash = ();
#   open(KNOWN, $file);
#   while(<KNOWN>) {
#     chomp;
#     #my($gb, $desc, @rest) = split("\t", $_);
#     my(@info) = split("\t", $_);
#     my $gb = $info[0];

#     # debug
#     #next unless $gb eq "BA000030";
#     #warn "gb eq $gb, line = '$_'\n";

#     my $cluster = $info[4]; # new
#     $cluster =~ /Cluster (\d+)/;
#     my $clusternum = $1;

# # New Oct 2013
#     if($_ =~ /$gene_cluster_phrase/i) { 
# 	#print "FOUND phrase for $gb\n"; # debug
# 	# get_gene_cluster_phrase($info[1]) || 
# 	my $phrase = get_gene_cluster_phrase($info[19]) || get_gene_cluster_phrase($info[20]) || "Warning, gene cluster phrase found but description not found";
# 	#print "phrase = $phrase\n"; # debug
# 	$hash{$gb}{$clusternum} = $phrase; 
#     }

###### Old before paper revisions:  was too specific, missed some known ones?  This led to < 172
#    my $desc = $info[1];
#    my $synline = $info[18] || "";
#    my $simline = $info[19] || "";
#    my $phrase = "";
#    if ($phrase = get_gene_cluster_phrase($desc)) {$hash{$gb} = $phrase;}
#    elsif ($phrase = get_gene_cluster_phrase($synline)) {$hash{$gb} = $phrase;}
#    elsif ($phrase = get_gene_cluster_phrase($simline)) {$hash{$gb} = $phrase;}
#    # Entire synonym - too much

#     else {$hash{$gb}{$clusternum} = 0;}
#   }
#   close(KNOWN);
#   return (%hash);
# }
# sub get_gene_cluster_phrase {
#   my $descline = shift;

#   # debug
#   #print "descline = $descline\n";
 
#   if($descline !~ /;/){ # only one genbank
#     if ($descline =~ /.*bp (.*$gene_cluster_phrase.*)/i) { return $1; }
#     #else {return $descline;}
#   }
#   #if($descline =~ /;/) { # multiple genbanks here
#   else {
#     #print "HERE found multiple\n";
#     my(@descs) = split(';', $descline);
#     foreach my $desc(@descs) { 
#       #print "  HERE desc = $desc\n";
#       if ($desc=~ /.*bp (.*$gene_cluster_phrase.*)/i) { 
# 	#print "    HERE FOUND, returning '$1'\n";
# 	return $1; 
#       }
#     }
#   }
#   return 0;
# }
