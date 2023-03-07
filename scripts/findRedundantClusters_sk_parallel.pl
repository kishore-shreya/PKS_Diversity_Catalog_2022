#!/usr/bin/perl -w
use strict;

die "Usage:  cluster_file antismash_dir start end\n" unless @ARGV >= 4;

my $infile = shift;
my $antismash_dir = shift;
my $line_num_start = shift; #0 is the first line
my $line_num_end = shift; # iterate until this line (inclusive)
my $end_of_line = shift; #the number of lines - 1
my $dest_file = shift;

warn "Building GB sequences...\n";
my %gb_clusternum_seq = build_gb_seq_dict($infile);
my @lines_range = ($line_num_start .. $line_num_end);

warn "Iterating to compare for identical sequences...\n";

open (my $fh, '>', $dest_file);

foreach my $line_num (@lines_range){
    my $cur_accession = $gb_clusternum_seq{$line_num}{"acc"};
    my $cur_clusternum = $gb_clusternum_seq{$line_num}{"num"};
    my $cur_coords = $gb_clusternum_seq{$line_num}{"coords"};
    my $cur_seq = $gb_clusternum_seq{$line_num}{"seq"};
    my ($cur_start, $cur_stop) = split('-', $cur_coords);
    my $cur_length = $cur_stop - $cur_start;
    #warn "Printing the current line number $line_num and the accession $cur_accession";
    my @second_range = ($line_num .. $end_of_line);
    foreach my $other_line_num (@second_range){
        my $other_accession = $gb_clusternum_seq{$other_line_num}{"acc"};
        my $other_clusternum = $gb_clusternum_seq{$other_line_num}{"num"};
        my $other_coords = $gb_clusternum_seq{$other_line_num}{"coords"};
        my $other_seq = $gb_clusternum_seq{$other_line_num}{"seq"};
        my ($other_start, $other_stop) = split('-', $other_coords);
        my $other_length = $other_stop - $other_start;
	#warn "going outside  $cur_accession $other_accession $line_num $other_line_num";
	#warn "printing $other_coords";
        #warn "Printing the other line number $other_line_num and the accession $other_accession";
        
        if($cur_seq eq $other_seq) { 
	  #warn "sequence are equal";
	  #print "hello in middle";
          print $fh "$cur_accession\t$cur_clusternum\t$cur_coords\t$other_accession\t$other_clusternum\t$other_coords\t$cur_length\t$other_length\n";
        }
	   elsif($cur_seq =~ /$other_seq/ || $other_seq =~ /$cur_seq/) { ############################################## one is a substr of the other
          print $fh "$cur_accession\t$cur_clusternum\t$cur_coords\t$other_accession\t$other_clusternum\t$other_coords\t$cur_length\t$other_length\tsubstring match\n";
        }  
    }
}


close $fh;

exit;

sub build_gb_seq_dict {
  my $infile = shift;
  my %gb_seq = ();
  my $line_num = 0;
  open(GBSEQ, $infile);
  while(<GBSEQ>) {
    chomp;
    my($accession, $desc, $clusternum, $coords, $seq) = split("\t", $_);
    if(!$seq) {
      die "Error, no seq for $accession\n";
    } 
    #warn "$clusternum";
    $gb_seq{$line_num} = {"acc"=>$accession, "num"=>$clusternum, "coords"=>$coords, "seq"=>$seq};
    $line_num++;
  }
  return %gb_seq;
}
