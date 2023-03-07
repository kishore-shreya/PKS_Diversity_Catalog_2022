#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
#jake hsu edit
#commented out novel_domain_file input 
#couldn't figure out how to generate it
#and it doesn't appear to get used to generate the output in this script
die "Usage: $0 cluster_file antismash_dir genbank_dir date_file sequence_file novel_domain_file\n" unless @ARGV >= 5;
my $cluster_file = shift;
my $antismash_dir = shift;
my $genbank_dir = shift;
my $date_file = shift; # not using
my $sequence_file = shift;
#my $novel_domain_file = shift;
my $url = "http://proteus.stanford.edu/~maureenh/antismash2.6dbs.min3.pksnrps";
my $html_file = "index.html";
my $js_filename = "geneclusters.js";
my @interesting_antismash_domains = qw/oMT nMT cMT Aminotran Hetero/;

my $sequence_error = "Error in sequence file";
#my $sequence_error = "NA";

warn "Building info...\n";
#my %accession_year = build_accession_year($date_file);
#my %accession_clusternum_domain = build_novel_domains($novel_domain_file);
my %accession_clusternum_sequence = build_accession_clusternum_sequence($sequence_file);

# Printing file...\n";
# Print header
print "NCBI sequence accession";
print "\tSequence description";
print "\tDate in NCBI";
print "\tSpecies";
print "\tCluster number within this sequence";
print "\tLink to graphical antismash view (click on cluster number on page)";
print "\tCluster coord-inates";
#print "\tNCBI sequence length (bp)";
print "\tCluster length (bp)";
print "\tCluster GC content (%)";
#print "\tNum orfs in cluster";
#print "\tNum PKS orfs in cluster";
print "\tNum KSs";
print "\tNum ATs";
print "\tNum ACPs";
print "\tNum KRs";
print "\tNum DHs";
print "\tNum ERs";
print "\tNum Cs";
print "\tNum As";
print "\tInteresting Antismash-predicted domains"; # skip this because of bug
#print "\tUnusual inter-domain gaps (from blast NR)";
print "\tCluster type";
print "\tClusters with identical sequence";
print "\n";



warn "Reading from $cluster_file\n";
open(IN, $cluster_file);
while(<IN>) {

  chomp;
  my($accession, $desc, $year, $clusternum, $coords, $length, $synonyms) = split("\t", $_);
  $clusternum =~ m/\"r1c(.*)\"/;
  my $size = length $1;
  my $result_string;
  if ($size == 1){$result_string = "00$1";}
  if ($size == 2){$result_string = "0$1";}
  if ($size == 3){$result_string = "$1";}

  my $antismash_accession = $accession;
  if($accession =~ /(.*)\.\d+/) {
    $antismash_accession = $1;
  }
  my $local_antismash_image_file;
  if (-e "$antismash_dir/$antismash_accession/$antismash_accession.region$result_string.gbk")
  {	
	$local_antismash_image_file = "$antismash_dir/$antismash_accession/$antismash_accession.region$result_string.gbk";
  }
  elsif (-e "$antismash_dir/$antismash_accession/$antismash_accession.1.region$result_string.gbk")
  {
	$local_antismash_image_file = "$antismash_dir/$antismash_accession/$antismash_accession.1.region$result_string.gbk";
  }
  my $hyperlink = "$url/$accession/$html_file#cluster-$clusternum";
  my $genecluster_file = "$antismash_dir/$antismash_accession/geneclusters.txt";
  my $nrpspks_prots_file = "$antismash_dir/$antismash_accession/nrpspks_proteins.fasta";

  # for num domains, should be using embl_lines?
  #eg antismash-1.2.2.PKSonly.min4.4dbs/batch_output/AEFS01000009/embl_lines.txt
  # grep PKS_KS antismash-1.2.2.PKSonly.min4.4dbs/batch_output/AEFS01000009/embl_lines.txt
  my $num_kss = get_num_domains($local_antismash_image_file, "PKS_KS"); #HERE
  my $num_ats = get_num_domains($local_antismash_image_file, "PKS_AT");
  my $num_acps = get_num_domains($local_antismash_image_file, "ACP");
  my $num_krs = get_num_domains($local_antismash_image_file, "PKS_KR");
  my $num_dhs = get_num_domains($local_antismash_image_file, "PKS_DH");
  my $num_ers = get_num_domains($local_antismash_image_file, "PKS_ER");
  my $num_cs = get_num_Cdomains($local_antismash_image_file, "Domain: Condensation_");
  my $num_as = get_num_domains($local_antismash_image_file, "AMP-binding");

  my @type = get_cluster_type_new($local_antismash_image_file); 
  my($start, $stop) = split('-', $coords);
  my $length2 = $stop - $start;
  # my (@orfs) = get_orfs($genecluster_file, $clusternum);
  # my $num_orfs = @orfs;

  #my(@nrpspks_prots_total) = get_nrpspks_prots($nrpspks_prots_file);
  #my @nrpspks_orfs  = overlap(\@nrpspks_prots_total, \@orfs);
  #my $num_nrpspks_orfs = @nrpspks_orfs;

  my $species = get_species($desc);
  #my $year = $accession_year{$accession};
  #my ($gc, $total_genbank_len) = get_gc_content_length($accession, $coords);
  my ($gc) = get_gc_content($accession, $clusternum);
  unless ($gc =~ /Error/) {
  $gc = $gc *100;
  $gc = sprintf("%.2f", $gc);
  }

  my $antismash_domains = ""; #get_antismash_interesting_domains($accession); # TODO pass clusternum and fix this if want to include.  Right now has a bug where it pulls all domains from the genbank id
  #my $blasthits = $accession_clusternum_domain{$accession}{"c$clusternum"};
  #my (@novel_domains) = sort keys %$blasthits;
  #my $novel_domains = join(';', @novel_domains);
  

  print "$accession";
  print "\t$desc";
  print "\t$year";
  print "\t$species";
  print "\tCluster $clusternum";
  print "\t$hyperlink";
  print "\t$coords";
  print "\t$length";
  print "\t$gc";
  #print "\t$num_orfs";
  #print "\t$num_nrpspks_orfs";
  print "\t$num_kss";
  print "\t$num_ats";
  print "\t$num_acps";
  print "\t$num_krs";
  print "\t$num_dhs";
  print "\t$num_ers";
  print "\t$num_cs";
  print "\t$num_as";
  print "\t$antismash_domains";
  #print "\t$novel_domains";
  print "\t@type";
  print "\t$synonyms";
  print "\n";


}
exit;

#sub build_novel_domains {
#  my $file = shift;
# my %accession_clusternum_domain = ();
#  open(NOVEL, $file) || warn "Warning no novel domain file $file\n";
#  while(<NOVEL>) {
#    my($accession, $clusternum, $interestingWord, $blasthit, $blastfile) =split("\t", $_);
#    $accession_clusternum_domain{$accession}{$clusternum}{$blasthit} = 1;
    #print "$accession\t$clusternum\t$blasthit\n";
#  }
  #exit;
  #return %accession_clusternum_domain;
#}

sub get_antismash_interesting_domains {
  my $accession = shift;
  my $file = "$antismash_dir/$accession/$js_filename";
  my $antismash_domains = "";
  #open(ANTI, $file) || die "Can't open $file\n";;
  #while(<ANTI>) {
  #  next unless $_ =~ /misc_feature/;
  #  my($a, $b, $domain, @rest) = split("\t", $_);
    foreach my $interesting_domain(@interesting_antismash_domains) {
  #    if($_ =~ /$interesting_domain/) {
      my $result = `grep $interesting_domain $file`;
      if ($result) {
        #chomp $result;
	#$result =~ s/\s//g;
        $antismash_domains .= " $interesting_domain";
      }
    }
  #}
  #close(ANTI);
  return $antismash_domains;
}
sub get_gc_content_length_old {
  my $gb = shift;
  my $coords = shift;
  #my $file = "$antismash_dir/$gb/$gb.fasta";
  my $file = "$genbank_dir/$gb.gb";
  my $format = "genbank";

  unless (-e $file) { die "Error, no file $file\n";}
  my $inseq = Bio::SeqIO->new(
                            -file   => "<$file",
                            -format => $format,
                            );

  my $seqobj = $inseq->next_seq;
  my $seq = $seqobj->seq;
  my $total_length = $seqobj->length; ################################
  my($start, $stop) = split('-', $coords);
  my $len_start_stop = $stop - $start;
  my $cluster = substr($seq, $start, $len_start_stop);
  #print "cluster= $cluster\n";
  my $len = length $cluster;
  my $len_seq = length $seq;
  my $gcount = ($cluster =~ tr/G//);
  my $ccount = ($cluster =~ tr/C//);

  my $gc_num = $gcount+$ccount;
  my $gc = $gc_num / $len;
  return ($gc, $total_length);
}



sub get_gc_content {
  my $accession = shift;
  my $clusternum= shift;
  my $seq = $accession_clusternum_sequence{$accession}{$clusternum};

  my $len_seq = length $seq;
  my $gcount = ($seq =~ tr/G//);
  my $ccount = ($seq =~ tr/C//);

  my $gc_num = $gcount+$ccount;
  my $gc = $gc_num / $len_seq;
  return $gc;
}

sub get_species {
  my $desc = shift;
  my(@words) = split(' ', $desc);
  my $genus = $words[0];
  my $species = $words[1];
  if($species =~ /sp\./) {
    $species .= " $words[2]";
    if($words[2] =~ /ATCC/) {
      $species .= " $words[3]";
    }
  }
  return "$genus $species";
}

sub overlap {
  my $orfs1 = shift;
  my $orfs2 = shift;
  my @orfs = ();
  foreach my $orf1(@$orfs1) {
    foreach my $orf2(@$orfs2) {
      if($orf1 eq $orf2) {
        push(@orfs, $orf1);
      }
    }
  }
  return (@orfs);
}


sub get_nrpspks_prots {
  my $file = shift;
  open(PROT, $file) || die "Can't open file $file\n";
  my @orfs = ();
  while(<PROT>) {
    chomp;
    if($_ =~ /\>(.*)/) {
      push(@orfs, $1);
    }
  }
  close(PROT);
  return (@orfs);
}

sub get_orfs{
  my $file = shift;
  my $clusternum = shift;
  my $cluster = "c$clusternum";
  my $num_orfs = 0;
  warn "$file with issues is";
  open(TRANS, $file);
  while(<TRANS>) {
    chomp;
    my ($gb, $desc, $cluster2, $type, $orfs1, $orfs2) = split("\t", $_);
    if ($cluster eq $cluster2) {
      my(@orfs) = split(';', $orfs1);
      return (@orfs);
    }
  }
  close (TRANS);
}


sub get_cluster_type_new {
  my $file = shift;
  my $grep_query = "\/product";
  my $results = `grep $grep_query $file`;
  my @results = split("\n", $results);
  my %counts;
  foreach my $result (@results) {
    $result=~ m/$grep_query.*\"(.*)\"/;
    $counts{$1} = 1;
  }
  return keys %counts;
}

sub get_cluster_type{
  my $file = shift;
  my $clusternum = shift;
  #my $cluster = "c$clusternum";
  my $type = "";
  open(TRANS, $file);
  my $file_clusternum = 1;
  while(<TRANS>) {
    chomp;
    my ($gb, $desc, $type, @genes) = split("\t", $_);
    if ($clusternum eq $file_clusternum) {
      return $type;
    }
    $file_clusternum++;
  }
  close (TRANS);
}

sub get_num_domains {
  my $file = shift;
  my $domain = shift;
  my $grep_query = "aSDomain.*\"$domain\"";
  my $results = `grep $grep_query $file`;
  my @results = split("\n", $results);
  my $counts=0;
  foreach my $result (@results) {
    $counts++;
  }
  return $counts;
}

sub get_num_Cdomains {
  my $file = shift;
  my $domain = shift;
  my $grep_query = "Domain:.*Condensation_";
  my $results = `grep $grep_query $file`;
  my @results = split("\n", $results);
  my $counts=0;
  foreach my $result (@results) {
    $counts++;
  }
  return $counts;
}


sub get_num_KSs_old {
  my $file = shift;
  my $results = `grep 'KS</text>' $file`;
  my @results = split("\n", $results);
  my $num_kss = @results;
  if($num_kss ==0) {
    $results = `grep KS_domain $file`;
    @results = split("\n", $results);
    $num_kss = @results;
  }
  return $num_kss;
}


sub get_num_DHs_old {
  my $file = shift;
  my $results = `grep 'DH</text>' $file`;
  my @results = split("\n", $results);
  my $num_kss = @results;
  if($num_kss ==0) {
    $results = `grep DH_domain $file`;
    @results = split("\n", $results);
    $num_kss = @results;
  }
  return $num_kss;
}

sub get_num_ERs_old {
  my $file = shift;
  my $results = `grep 'ER</text>' $file`;
  my @results = split("\n", $results);
  my $num_kss = @results;
  if($num_kss ==0) {
    $results = `grep ER_domain $file`;
    @results = split("\n", $results);
    $num_kss = @results;
  }
  return $num_kss;
}


sub get_num_ATs_old {
  my $file = shift;
  my $results = `grep 'AT</text>' $file`;
  my @results = split("\n", $results);
  my $num_kss = @results;
  if($num_kss ==0) {
    $results = `grep AT_domain $file`;
    @results = split("\n", $results);
    $num_kss = @results;
  }
  return $num_kss;
}


sub build_accession_clusternum_sequence {
  my $file = shift;
  open(SEQ, $file);
  my %hash = ();
  while(<SEQ>) {
    my($accession, $desc, $clusternum, $coords, $seq) = split("\t", $_);
    $hash{$accession}{$clusternum} = $seq;
  }
  return (%hash);
}

