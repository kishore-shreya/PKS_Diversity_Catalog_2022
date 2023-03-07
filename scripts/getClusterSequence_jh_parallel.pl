#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

die "Usage: $0 genbank_dir cluster_file antismash_dir start end\n" unless @ARGV >= 3;
my $genbank_dir = shift;
my $cluster_file = shift;
my $antismash_dir = shift;
my $pair_count_start = shift;
my $pair_count_end = shift || 0;
my $sequence_error = "Error in sequence file";
my $pair_counter = 0;

open(IN, $cluster_file);
while(<IN>) {

  chomp;
  my($accession, $desc, $clusternum, $coords) = split("\t", $_);
  if($pair_counter < $pair_count_start) {$pair_counter++; next;}############################## parallelize
        if($pair_count_end){
          if($pair_counter >= $pair_count_end) {exit;}
        }
  my $antismash_accession = $accession;
  if($accession =~ /(.*)\.\d+/) {
    $antismash_accession = $1;
  }
  my $local_antismash_image_file = "$antismash_dir/$antismash_accession/svg/genecluster$clusternum.svg";

  my($start, $stop) = split('-', $coords);
  my $length = $stop - $start;

  my ($gene_cluster_seq) = get_cluster_seq($accession, $coords);
  print "$accession\t$desc\t$clusternum\t$coords\t$gene_cluster_seq\n";
  
  # calculate gc content
  #unless ($gene_cluster_seq=~ /Error/) {
  #$gc = $gc *100;
  #$gc = sprintf("%.2f", $gc);
  #}
  $pair_counter++;

}
exit;

sub get_cluster_seq{
  my $gb = shift;
  my $coords = shift;
  my $file = "$genbank_dir/$gb.gb";
  my $format = "genbank";

  unless (-e $file) { die "Error, no file $file\n";}
  my $inseq = Bio::SeqIO->new(
                            -file   => "$file",
                            -format => $format,
                            );

  my $seqobj = $inseq->next_seq;
  my $seq = $seqobj->seq;
  my $total_length = $seqobj->length; 
  my($start, $stop) = split('-', $coords);
  my $len_start_stop = $stop - $start;
  my $cluster = substr($seq, $start, $len_start_stop);
  return $cluster;
}

#  my $len = length $cluster;
#  my $len_seq = length $seq;
#  my $gcount = ($cluster =~ tr/G//);
#  my $ccount = ($cluster =~ tr/C//);
#  my $gc_num = $gcount+$ccount;
#  my $gc = $gc_num / $len;
#  return ($gc, $total_length);
#}




