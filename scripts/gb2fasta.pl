#!/usr/bin/perl -w
# 
#
use strict;
use Bio::SeqIO;

die "Usage $0 genbank fasta\n" unless @ARGV >= 2;

my $genbank_file = shift;
my $fasta_file = shift;

print_fasta($genbank_file, $fasta_file);

sub print_fasta {
  my $genbank_file = shift;
  my $fasta_file = shift;

  my $in  = Bio::SeqIO->new(-file => $genbank_file ,
                           -format => 'genbank');
  my $out = Bio::SeqIO->new(-file => ">$fasta_file" ,
                           -format => 'Fasta');

  my $seq = $in->next_seq();
  $out->write_seq($seq);
}

