#!/usr/bin/perl -w
# 
# Choose one synonym, skip the rest.
# Choose the first for now.
# Print the synonymous accession #s at the end.
#

#Jake Hsu edit
#comment out lines appending synonyms descriptions (too confusing)
#comment out lines getting new species (too confusing)


use strict;

die "Usage $0 tableFile similarFile accessionFile dateFile\n" unless @ARGV >= 4;

my $table_file = shift;
my $similar_file = shift;
my $accession_desc_file = shift;
my $date_file = shift;

my %gb_synonyms = build_gb_synonyms($similar_file); # key = accession_clusternum1, key2 = accession_clusternum2

my %accession_desc = build_accession_desc_no_version($accession_desc_file);
my %accession_year = build_accession_year($date_file);

# Oct 2013
my %accession_clusternum_length = build_accession_clusternum_length($table_file);

my %skip = ();



open(IN, $table_file);
my $header = <IN>;
chomp $header;
print "$header\tClusters with > 90% similarity score (chained through transitivity)\n";

while(<IN>) {
  chomp;
  my(@info) = split("\t", $_);
  # Need accession and clusternum
  my $accession = $info[0];
  my $year = $info[2];
  my $clusternum = $info[4];
  my $length = $info[7];
  my $accession_clusternum = "$accession\t$clusternum";
  my $desc = $accession_desc{$accession} || warn "Can't find desc for $accession\n";

  if($skip{$accession_clusternum}) {next;} #################################################### ?
  #my $longest_length = $length; # if use the check below


  # Get synonyms for this gb (from synonym file)
  my $synonyms = $gb_synonyms{$accession_clusternum};
  #my $year = $accession_year{$accession};
  my $synonym_str = "$accession Date=$year $clusternum $length"."bp"." $desc";

  # Iterate synonyms and choose first, make a note of it
  foreach my $synonym (keys %$synonyms) {
    my($syn_accession, $syn_clusternum) = split("\t", $synonym);
    next if($accession eq $syn_accession);
    $skip{$synonym} = 1; #################################################### ?
    my $syn_desc = $accession_desc{$syn_accession} || warn "Can't find desc for $syn_accession\n";
    my $syn_year = $accession_year{$syn_accession} || warn "Can't find year for $syn_accession\n";
    my $syn_length = $accession_clusternum_length{$syn_accession}{$syn_clusternum} || warn "Can't find length for $syn_accession\n";

# TODO: add this check/edit to change the accession number?
#    if($syn_length > $longest_length) {
#      $longest_length = $syn_length;
#      $print_accession = $syn_accession; # if this synonym is longest, print this one instead
#      $print_clusternum = $syn_clusternum;
#      $print_desc = $syn_desc;
#      $print_coords = $syn_coords;
#    }


    # Is this year the oldest?  If so, reset the year
    if($syn_year < $info[2]) {
      $info[2] = $syn_year;
    #  $info[1] = "$info[1]; $syn_desc";
    }

    #$synonym_str .= "; $syn_accession "."cluster$syn_clusternum $syn_length"."bp"." $syn_desc"; # Add the info for this synonym to the string
    #$synonym_str .= "; $syn_accession Date=$syn_year "."$syn_clusternum $syn_desc"; # Add the info for this synonym to the string
    $synonym_str .= "; $syn_accession Date=$syn_year "."$syn_clusternum $syn_length"."bp"." $syn_desc"; # Oct 2013

  }

  push(@info, $synonym_str);
  # Check that species has a species
  #if($info[3] =~ /Sequence \d/) {
  #   $info[3] = get_species_from_info(@info) || $info[3];
  #}
  my $line_new = join("\t", @info);
  #print "$line_new\t$synonym_str\n";
  print "$line_new\n";
  #print "$_\t$synonym_str\n";
}

exit;

sub get_species_from_info{
  my(@info) = @_;

  for(my $i=19; $i<=21; $i++) {
    my(@descs) = split(';', $info[$i]);
    #print "i = $i, descs = '@descs'\n";
    foreach my $desc(@descs) {
      my $species = get_species_from_desc($desc);
      if($species =~ /Sequence \d/ || !$species) {next;}
      return $species;
    }
  }
  return 0;
}

sub get_species_from_desc {
  my $desc = shift;
  $desc =~ /bp (.+? .+? )/;
  return $1 || 0;
}

sub build_gb_synonyms {
  my $infile = shift;
  my %synonyms = ();
  open(GB, $infile);
  while(<GB>) {
    chomp;
    my($accession1, $clusternum1, $accession2, $clusternum2, @rest) = split("\t", $_);
    my $accession_clusternum1 = "$accession1\t$clusternum1";
    my $accession_clusternum2 = "$accession2\t$clusternum2";

    # skip self-self ones?
    $synonyms{$accession_clusternum1}{$accession_clusternum2} = 1;
    $synonyms{$accession_clusternum2}{$accession_clusternum1} = 1;
  }
  close(GB);
  return %synonyms;
}

sub build_accession_desc_no_version{
  my($file) = shift;
  open(TWOCOL, $file) || die "Error can't open $file\n";
  my %hash;
  while(<TWOCOL>) {
    chomp;
    # debug
    #if($_ =~ /BSU11039/) {print "here line = $_\n";}
    my($first, $second, @rest) = split("\t", $_);
    my $accession = $first;
    if($first =~ /(.*)\.\d/) {
      $accession = $1;
    }
    #if($_ =~ /BSU11039/) {print "here accession = '$accession'\n";}
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


sub build_accession_clusternum_length {
  my $file = shift;
  my %hash = ();
  open(LEN, $file);
  while(<LEN>) {
    chomp;
    my(@info) = split("\t", $_);
    $hash{$info[0]}{$info[4]} = $info[7];
  }
  return (%hash);
}


sub print_hash{
  my(%hash) = @_;
  foreach my $key(sort keys %hash) {
    print "'$key'\t'$hash{$key}'\n";
  }
}


