#!/usr/bin/perl -w
# 
# Choose one synonym, skip the rest.
# Choose the longest for now.
# Print the synonymous accession #s at the end.
#

use strict;

die "Usage:  $0 infile synonym_file antismash_dir accession_desc_file\n" unless @ARGV >= 4;
my $infile = shift;
my $synonym_file = shift;
my $antismash_dir = shift;
my $accession_desc_file = shift;
my $date_file = shift;

my %gb_synonyms = build_gb_synonyms($synonym_file); # key = accession_clusternum1, key2 = accession_clusternum2
my %accession_desc = build_accession_desc_no_version($accession_desc_file);
my %accession_year = build_accession_year($date_file);
my %skip = ();
my %accession_clusternum_coords = build_accession_clusternum_coords($infile);

#my %accession_print_accession = (); # TODO
#my %accession_synonym_str = (); # TODO
my %seen = ();

# First choose which synonym to print
# Then print

open(IN, $infile);
while(<IN>) {
  chomp;
  my($accession, $desc, $clusternum, $coords) = split("\t", $_);

  $desc =~ s/;/,/g; # Oct 2013

  my $accession_clusternum = "$accession\t$clusternum";
  #my $coords = $accession_clusternum_coords{$accession}{$clusternum};
  my($start, $stop) = split('-', $coords);
  my $length = $stop - $start;

  if($skip{$accession_clusternum}) {next;} ##

  my $longest_length = $length;
  my $print_accession = $accession; # At first, select this accession as the synonym to print
  my $print_desc = $desc;
  my $print_clusternum = $clusternum;
  my $print_coords = $coords;
  my $year = $accession_year{$accession};

  # Get synonyms for this gb (from synonym file)
  my $synonyms = $gb_synonyms{$accession_clusternum};
  my $synonym_str = "$accession Date=$year "."cluster$clusternum $length"."bp"." $desc";
 
  # Iterate synonyms and choose longest, make a note of it
  foreach my $synonym (keys %$synonyms) {
    my($syn_accession, $syn_clusternum) = split("\t", $synonym);
    my $syn_desc = $accession_desc{$syn_accession};
    $syn_desc =~ s/;/,/g; # Oct 2013

    my $syn_year = $accession_year{$syn_accession};
    my $syn_coords = $accession_clusternum_coords{$syn_accession}{$syn_clusternum};
    my($syn_start, $syn_stop) = split('-', $syn_coords);
    my $syn_length = $syn_stop - $syn_start;
    if($syn_length > $longest_length) {
      $longest_length = $syn_length;
      $print_accession = $syn_accession; # if this synonym is longest, print this one instead
      $print_clusternum = $syn_clusternum;
      $print_desc = $syn_desc;
      $print_coords = $syn_coords;
    }
    if($syn_year < $year) {
      $year = $syn_year;
    }

    unless($accession eq $syn_accession) {
      $synonym_str .= "; $syn_accession Date=$syn_year "."cluster$syn_clusternum $syn_length"."bp"." $syn_desc"; # Add the info for this synonym to the string
    }

    $skip{$synonym} = 1; ##
  }
  if($seen{$print_accession}{$print_clusternum}) {next;} ##
  print "$print_accession\t$print_desc\t$year\t$print_clusternum\t$print_coords\t$longest_length\t$synonym_str\n";
  $seen{$print_accession}{$print_clusternum} = 1;
}



exit;




sub build_gb_synonyms {
  my $infile = shift;
  my %synonyms = ();
  open(GB, $infile);
  while(<GB>) {
    chomp;
    my($accession1, $clusternum1, $coords1, $accession2, $clusternum2, $coords2, $length1, $length2, @rest) = split("\t", $_);
    my $accession_clusternum1 = "$accession1\t$clusternum1";
    my $accession_clusternum2 = "$accession2\t$clusternum2";

    # skip self-self ones?
    $synonyms{$accession_clusternum1}{$accession_clusternum2} = 1;
    $synonyms{$accession_clusternum2}{$accession_clusternum1} = 1;
  }
  return %synonyms;
}

sub build_accession_clusternum_coords {
  my $file = shift;
  my %hash = ();
  open(LEN, $file);
  while(<LEN>) {
    chomp;
    my($accession, $desc, $clusternum, $coords) = split("\t", $_);
    $hash{$accession}{$clusternum} = $coords;
  }
  return (%hash);
}

sub build_accession_desc_no_version{
  my($file) = shift;
  open(TWOCOL, $file) || die "Error can't open $file\n";
  my %hash;
  while(<TWOCOL>) {
    chomp;
    my($first, $second, @rest) = split("\t", $_);
    $first =~ /(.*)\.\d/;
    my $accession = $1 || $first;
    $hash{$accession} = $second;
  }
  close(TWOCOL);
  return %hash;
}


sub build_accession_year {
  my $file = shift;
  my %hash = ();
  my %accession_date = build_hash_from_two_cols($file);
  foreach my $accession(keys %accession_date) {
    my $date = $accession_date{$accession};
    $date =~ /(.*)\-(.*)\-(.*)/;
    my $year = $3 || "Error in year";
    $hash{$accession} = $year;
  }
  return (%hash);
}

sub build_hash_from_two_cols {
  my($file) = @_;
  open(TWOCOL, $file) || die "Can't open $file\n";
  my %hash;
  while(<TWOCOL>) {
    chomp;
    my($first, $second, @rest) = split("\t", $_);
    $hash{$first} = $second;
  }
  return %hash;
}


