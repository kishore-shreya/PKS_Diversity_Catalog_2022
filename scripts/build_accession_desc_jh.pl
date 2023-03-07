#!/usr/bin/perl -w
use strict;

die "Usage: $0 infile\n" unless @ARGV >= 1;

my $infile = shift;


  open(IN, $infile) || die "Can't open file $infile\n";
  while (<IN>) {
    chomp;
    my($genbank, $desc, @rest) = split("\t", $_);
    #unless($genbank =~ /.*\|(.*)\|(.*)/) {
      #die "Error, not two pipes in $genbank\n";
    #}
    #$genbank =~ s/ref//g;
    #$genbank =~ s/dbj//g;
    #$genbank =~ s/emb//g;
    #$genbank =~ s/gb//g;
    #$genbank =~ s/tpe//g;
    #$genbank =~ s/\|//g;
    #print "$genbank\t$desc\n";
    print "$genbank\t$desc\n";

    # Some lines have a second accession... print these
    if($2) { print "$2\t$desc\n";}
  }
  close(IN);




