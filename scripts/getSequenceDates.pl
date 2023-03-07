#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $folder = shift;
my $format = "genbank";
 
my @files = `ls $folder/*gb`;
chomp(@files);
foreach my $file(@files){
  #print "file = $file\n";

  my $inseq = Bio::SeqIO->new(
                            -file   => "$file",
                            -format => $format,
                            );

  my $seq = $inseq->next_seq;
  my $accession = $seq->accession_number;
  

  $file =~ /$folder\/(.*)\.gb/;
  my $filename = $1 || warn "warning, no filename for $file\n";
  #warn $seq->accession_number, "\n";
  unless ($accession eq $filename) {
    warn "warning, filename $filename ne accession $accession...\n";
  }

  print "$filename";
  my (@dates) =  $seq->get_dates;
  foreach my $date(@dates) {
    print "\t$date";
  }

  my $anno_collection = $seq->annotation; 
  my @annotations = $anno_collection->get_Annotations('reference');
  for my $value ( @annotations ) {
    my $hash_ref = $value->hash_tree;
    for my $key (keys %{$hash_ref}) {
      next unless $key eq 'location';
      my $location = $hash_ref->{$key};
      if($location =~ /(\d{2}-\D{3}-\d{4})/) { # Find date-formatted parts of this reference location
        print "\t$1";
      }
    }
  }
  print "\n";
}


