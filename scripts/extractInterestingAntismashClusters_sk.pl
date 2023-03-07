#!/usr/bin/perl -w
#shreya kishore edit
#modified to be able to read newer antiSMAHSH 6.1.1 output files
use strict;
#jake hsu edit
#This has been modified to find the start and end coordinates a little differently
# The start coordinates are the js file near the phrase "knowncluster" so look there
die "Usage: $0 accession_desc_file antismash_dir out_dir numKS\n" unless @ARGV >= 4;

my $accession_desc_file= shift;
my $antismash_batch_dir = shift;
my $output_folder = shift;
my $min_antismash_ks_num = shift;

# Outfile
my $outfile = "$output_folder/interestingClusters.txt";
my $no_results_outfile = "$output_folder/noAntismashResults.txt";
open(OUT, ">$outfile");
open(OUTNORESULT, ">$no_results_outfile");

# Accession description
my %accession_desc = build_hash_from_two_cols($accession_desc_file);
# Note that this should have 2826 lines
#print_hash(%accession_desc);
#exit;

my $debug = "SERERYAB";

foreach my $accession (sort keys %accession_desc) {
 
  # DEBUG
  #print "$accession\n"; # should lead to 2826 lines - yes. (2829 with extra at bottom)
  #next;
  #next unless $accession =~ /SER/; ########### debug


  #########################################################################
  # Antismash info
  # strip version
  $accession =~ /(.*)\.\d+/;
  my $antismash_accession = $1 || $accession;
  my $desc = $accession_desc{$accession};


  # DEBUG
  print "$antismash_accession\n";
  #next unless $antismash_accession eq $debug; ########### debug
  #if($antismash_accession eq $debug) {die "$debug\n";}


  #warn "$accession\n";
  my $genecluster_file = "$antismash_batch_dir/$antismash_accession/$antismash_accession.region001.gbk";
  my $genecluster_file_1 = "$antismash_batch_dir/$antismash_accession/$antismash_accession.1.region001.gbk";

  if (-e $genecluster_file) { ##################################################################?

    my $cmd = "ls $antismash_batch_dir/$antismash_accession/$antismash_accession.region*.gbk";
    ##print ".gbk files= $cmd\n";
    #warn "cmd = $cmd\n";

    my (@antismash_gbk_files) = `$cmd`;
  #if (@antismash_gbk_files) {
    chomp(@antismash_gbk_files);
    ##print "files = @antismash_gbk_files\n";

    my $num_files = @antismash_gbk_files;
    ##print "num files = $num_files\n";

    foreach my $x(@antismash_gbk_files) {
      #next; # WORKING!
      ##print "inside loop .gbk file = $x\n";
      my $num_kss = get_num_domains($x);
      ##print "number of KS = $num_kss\n";
      ##print "minimum KS = $min_antismash_ks_num\n";
      next unless $num_kss >= $min_antismash_ks_num;
      $x =~ /.*\/(.*\.gbk)/;
      my $filename = $1;
      $filename =~ /.*\.region0*([1-9]\d*)\.gbk/;
      my $clusternum = "\"r1c".$1."\"";
      ##print "clusternum = $clusternum";
      # get cluster coords from js file
      # eg antismash2.6dbs.min3.pksnrps/NC_000962/geneclusters.js 
      my $js_file = "$antismash_batch_dir/$antismash_accession/regions.js";
      #warn "js file = $js_file\n";
      my $coords = get_coords($js_file, $clusternum);
      print OUT "$antismash_accession\t$desc\t$clusternum\t$coords\n";
    }

    next;# get here ok.  Need to get this to stop working to figure out what's wrong.
  }
  elsif (-e $genecluster_file_1) {
      my $cmd = "ls $antismash_batch_dir/$antismash_accession/$antismash_accession.1.region*.gbk";
      ##print ".gbk files= $cmd\n";
      #warn "cmd = $cmd\n";

      my (@antismash_gbk_files) = `$cmd`;
    #if (@antismash_gbk_files) {
      chomp(@antismash_gbk_files);
      ##print "files = @antismash_gbk_files\n";

      my $num_files = @antismash_gbk_files;
      ##print "num files = $num_files\n";

      foreach my $x(@antismash_gbk_files) {
        #next; # WORKING!
        ##print "inside loop .gbk file = $x\n";
        my $num_kss = get_num_domains($x);
        ##print "number of KS = $num_kss\n";
        ##print "minimum KS = $min_antismash_ks_num\n";
        next unless $num_kss >= $min_antismash_ks_num;
        $x =~ /.*\/(.*\.gbk)/;
        my $filename = $1;
        $filename =~ /.*\.region0*([1-9]\d*)\.gbk/;
        my $clusternum = "\"r1c".$1."\"";
        ##print "clusternum = $clusternum";
        # get cluster coords from js file
        # eg antismash2.6dbs.min3.pksnrps/NC_000962/geneclusters.js
        my $js_file = "$antismash_batch_dir/$antismash_accession/regions.js";
        #warn "js file = $js_file\n";
        my $coords = get_coords($js_file, $clusternum);
        print OUT "$antismash_accession\t$desc\t$clusternum\t$coords\n";
      }

      next;# get here ok.  Need to get this to stop working to figure out what's wrong.
        
  }
  else { # no image files, no results
    print OUTNORESULT "No files for $accession in $antismash_batch_dir/$antismash_accession\n";
    next;
  }

}


print  "\nOutput to $outfile and $no_results_outfile\n\n";
exit;


# Subroutines

sub get_coords {
  my $js_file = shift;
  my $clusternum = shift;
  my $my_grep_line_start = "awk '/all_regions/,0' $js_file | grep '$clusternum' -m2 -A1 | tail -n1";
  my $my_grep_line_end = "awk '/all_regions/,0' $js_file | grep '$clusternum' -m2 -A2 | tail -n2";
  
  #awk 'f;/regexp/{f=1}' file
  ## print "\nmy_grep_line_start";
  ## print "\nmy_grep_line_end";
  #my $grep_search_actual_text = `$grep_start_search`;
  #print "\n$grep_search_actual_text";
  #my $grep_startline = "grep knowncluster --after-context=2 -m 1 $grep_search_actual_text";
  my $location_end_lines = `$my_grep_line_end`;
  my $location_start_lines =`$my_grep_line_start`;
  #print "\n$location_start_lines\n";
  #print "\n$location_end_lines\n";
  $location_start_lines =~ /start\": (.*),/;
  my $start = $1;
  $location_end_lines =~ /end\": (.*),/;
  my $end = $1;
  #warn "start = '$start', end = '$end'\n";
  my $coords = "$start-$end";
  #my $coords = 1-20;
  return $coords;
}
  


sub get_num_domains {
  my $file = shift;
  ## print "file=$file";
  my $grep_query = "aSDomain.*\"PKS_KS\"";
  ## print "grep_query=$grep_query";
  my $results = `grep $grep_query $file`;
  print "dollar result= $results";
  my @results = split("\n", $results);
  ## print "results= @results";
  my $counts=0;
  foreach my $result (@results) {
    $counts++;
  }
  return $counts;
}



sub build_hash_from_two_cols {
  my($file) = @_;
  open(TWOCOL, $file) || die "Can't open file $file\n";
  my %hash;
  while(<TWOCOL>) {
    chomp;
    my($first, $second, @rest) = split("\t", $_);
    $hash{$first} = $second;
  }
  return %hash;
}




sub print_hash{
  my(%hash) = @_;
  foreach my $key(sort keys %hash) {
    print "'$key'\t'$hash{$key}'\n";
  }
}

