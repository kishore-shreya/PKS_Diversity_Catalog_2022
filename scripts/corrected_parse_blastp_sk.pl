#!/usr/bin/perl -w

#sherya kishore edit to work with antiSMASH version 6 output files

# Copyright (C) 2019 Aleksandra Nivina

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

# Get % identity TOTAL in entire cluster
# = num matches in HSPs (non overlapping) / total length of cluster

# This script was extensively modified from the original script.
# Most importantly, it counts the hits in a way that avoids double-counting, which was the case for the original script.
# Also, it now works with 1-vs-all blastp or tblastx outputs, rather than 1-vs-1, as was the case for the original script.
# Now, you have to specify the range of such 1-vs-all output files (e.g. 0 1000; or 1000 2000), rather than the range of 1-vs-1 output files (e.g. 0 1000000; or 1000000 2000000)


use strict;

die "Usage: $0 infolder clusterfile outfile cluster_start_n cluster_end_n\n" unless @ARGV >= 5;

my $infolder = shift;
my $clusterFile = shift;
my $outfile = shift;
my $start = shift || 0;
my $end = shift || 0;


my (@files) = `ls $infolder`;
chomp(@files);
# 6987 files, of format ACCESSION.cluster00N.fasta.blastp.out

######### START OF MODIFIED PART ###########

my %query_lengths = build_query_length($clusterFile); # building a hash with keys: {accession}{clusternumber}, balues: cluster lengths

my %hits = ();
open(OUT, ">$outfile.$start-$end");

my @files_subset= @files[$start .. $end]; ############################## parallelize

#warn "Iterating files in $infolder\n";
my $counter = 0;
foreach my $file(@files_subset) { ############################## parallelize
  my $infile = "$infolder/$file";
  open(IN, $infile) || die "Error opening file $file\n";
  
  print "read file $infile\n";

  # get clusternum
  $file =~ /(.*)\.region(\d+)\.(.*)\.out/ || print "Filename for $file didn't have expected format\n"; # !!mod!!
  my $gb1 = $1;
  my $cluster1_prelim = $2;
  my $cluster1;
  my $first_letter = substr($cluster1_prelim, 0, 1);
  my $first_two_letters = substr($cluster1_prelim, 0, 2);

  if ($first_two_letters eq "00") {
  	$cluster1="r1c".substr($cluster1_prelim,2,1);
  }
  elsif ($first_letter eq "0") {
  	$cluster1="r1c".substr($cluster1_prelim,1,2);
  }
  else {
  	$cluster1="r1c".substr($cluster1_prelim,0,3);
  }

  # generating a variable that will store the current subject, and its hits
  my $prev_subj=""; 
  %hits = ();
  #warn "gb1 and cluster 1 are $gb1, $cluster1";
  my $query_length=$query_lengths{$gb1}{$cluster1};
  #print "query length is $query_length \n\n";
  
  while(<IN>) {
  	chomp;
  	next if $_ =~ /^\#/;
	next if $_ =~ /^$/;
  	my($query, $subject, $identity, $aln_length, $mismatches, $gapopens, $querystart, $queryend, $subjectstart, $subjectend, $evalue, $bitscore) = split("\t", $_);
	# print "read query $query, subject $subject\n";
	
  	# check if the subject is different from the one in previous line,
  	# in which case we save the score for the previous comparison and reset hits
  	if (("$subject" ne "$prev_subj") and ("$prev_subj" ne "")) {
  		# recording the score
		# print "preparing to record the score for the previous subject $prev_subj with the following info given to subprocess:\n";
	    my ($gb2,$cluster2)=split(/\./,$prev_subj);
	    $cluster2 =~ /\"(r1c.*)\"/;
	    my $size = keys %hits;
	    # print "query $gb1 $cluster1 of length $query_length\n";
	    # print "subject $gb2 $cluster2 and $size hit(s)\n";
	    record_score($gb1,$cluster1,$query_length,$gb2,$1,\%hits);
  		%hits = ();
  	}

  	my $queryspan = "$querystart\t$queryend";
	
  	# Swap start/end if necessary
  	if($querystart > $queryend) {
  		$queryspan = "$queryend\t$querystart";
  	}
	
	# print "query span $queryspan\n";
	
  	# Check whether an overlapping HSP has already been seen.   (Eg in tblastx, an HSP in a different frame.)
  	# The previous HSP had a better score (so blast returned it first).  So keep the % identity in the previous region the same.
  	# Need to get the % identity in this new overlapping region, excluding the previously-seen portion.
			
  	# checking if this query is a subset of any hit, includes or overlaps with any hit;
	# print "checking for overlap or include\n";
  	my @newspans=subset_overlap_or_include($queryspan,\%hits);
  	if (@newspans) {
  		foreach my $newspan (@newspans){
  			$hits{$newspan} = $identity;
  		}
  	}

  	# setting the current subject
  	$prev_subj=$subject;
  }
  
  # recording the score for the last subject
  if ($prev_subj eq ""){
  	print "here \n";
  }
  print "blastp $gb1 and $prev_subj \n";
  # print "preparing to record the score for the last subject with the following info given to subprocess:\n";
  my ($gb2,$cluster2)=split(/\./,$prev_subj);
  $cluster2 =~ /\"(r1c.*)\"/;
  my $size = keys %hits;
  # print "query $gb1 $cluster1 of length $query_length\n";
  # print "subject $gb2 $cluster2 and $size hit(s)\n";
  record_score($gb1,$cluster1,$query_length,$gb2,$1,\%hits);

  close(IN);
}

close(OUT);
exit;
  
######### END OF MODIFIED PART ###########

sub build_query_length {
  my $file = shift;
  my %hash = ();  # key = query, key = clusternum, val = length
  open(LEN, $file);
  while(<LEN>) {
    chomp;
    #my($gb, $desc, $clusternum, $coords) = split("\t", $_);
    my(@info) = split("\t", $_);
    my $gb = $info[0];
    #print "info = @info\ninfo0 = $gb";
    my $clusterstring = $info[4] || print "No value for clusterstring in line $_, file $file";
    $clusterstring =~ /Cluster\s\"(.+)\"/;
    my $clusternum = $1 || $clusterstring;
    my $len = $info[7];
    #my($start, $stop) = split("-", $coords);
    #my $len = $stop - $start;
    $hash{$gb}{$clusternum} = $len;
    #warn "$gb\t$clusternum\t$len\n";
  }
  close (LEN);
  return (%hash);
}

sub get_num_matches_bp {
  my($hits) = shift;
  my $num_matches_total = 0;
  foreach my $queryspan(keys %$hits) {
    my $identity = $$hits{$queryspan};
    my ($start, $end) = split("\t", $queryspan);
    my $len = $end - $start;
    if($len < 0) { ############################################ for hits on reverse strand
      die "Error:  length < 0 (length = $len, queryspan = $queryspan)\n\n";
      #$len = $len * -1;
    }
    my $num_matches = ($identity * $len ) / 100; # divide by 100 because % is given
    $num_matches_total += $num_matches;

    #print "queryspan $queryspan, num matches = $num_matches, identity = $identity, len = $len, num matches total = $num_matches_total\n";
  }
  return $num_matches_total;
}

######### START OF MODIFIED PART ###########

sub subset_overlap_or_include {
    my $queryspan = shift;
	my ($querystart, $queryend) = split("\t", $queryspan) || print "No query start and end for $queryspan\n";
	# print "Checking overlaps\n";
	my @extensions=$queryspan;
    my $hits = shift;
	foreach my $hit(keys %$hits){
		my ($hitstart, $hitend) = split("\t", $hit);
		# print "hit: $hit\n";
		my $index0=0;
		foreach my $extension(@extensions){
			# print "extension $extension:\n";
			my ($extensionstart, $extensionend) = split("\t", $extension);
			
			if ($extensionstart < $hitstart) {
				
				if ($extensionend > $hitstart) {
					
					# option #1: extension includes the hit
					if ($extensionend > $hitend){
						# print "includes the hit $hit\n";
						splice @extensions, $index0, 1;
						my $newspan1="$extensionstart\t$hitstart";
						push @extensions, $newspan1;
						my $newspan2="$hitend\t$extensionend";
						push @extensions, $newspan2;
					}
					
					# option #2: extension overlaps on the left
					else {
						# print "overlaps on the left with the hit $hit\n";
						splice @extensions, $index0, 1;
						my $newspan="$extensionstart\t$hitstart";
						push @extensions, $newspan;
					}				
				}
				
				# implicitly option #3: $extensionend <= $hitstart means that extension ends before the hitstart; don't do anything
				
			}
			
			elsif ($extensionstart < $hitend) {
				
				# option #4: extension is a subset of the hit
				if ($extensionend <= $hitend) {
					# print "is a subset of the hit $hit\n";
					splice @extensions, $index0, 1;
				}
				
				# option #5: extension overlaps on the right
				else {
					# print "is overlaps on the right with the hit $hit\n";
					splice @extensions, $index0, 1;
					my $newspan="$hitend\t$extensionend";
					push @extensions, $newspan;
				}
			}
				
			# implicitly option #6: $extensionstart>=$hitend means that extension starts after the hitend; don't do anything
						
			$index0++;
						
		}
	}
	
	# once we went over all extensions, check if some extensions are adjacent,
	# in which case we concatenate them into larger extensions
	@extensions=sort(@extensions);
	my $index1=0;
	my $m=@extensions-1;
	# print "Extentions for this query are: @extensions\n";
	# print "Checking if extensions can be concatenated\n";
	foreach my $extension1(@extensions){
		my ($extension1start,$extension1end)=split("\t", $extension1);
		for (my $index2=$index1+1; $index2 <= $m; $index2++) {
			my $extension2=$extensions[$index2];
			my ($extension2start,$extension2end)=split("\t", $extension2);
			if ($extension1end==$extension2start or $extension1end+1==$extension2start) {
				# print "Concatenated $extension1 and $extension2\n";
				$extensions[$index1]="$extension1start\t$extension2end";
				splice @extensions, $index2,1;
				$index2--;
				$m--;
				# print "new extensions: @extensions\n";
			}
		}
		$index1++;
	}
	
	return @extensions;
}

sub record_score {
	my $gb1 = shift;
	my $cluster1 = shift;
	my $query_length=shift;
	my $gb2=shift;
	my $cluster2=shift;
	my ($hits) = shift;
	# my $size=keys %hits;
	# print "preparing to record the score for the previous subject with the following info received by the subprocess:\n";
	#     print "subject $gb1 $cluster1 of length $query_length\n";
	#     print "query $gb2 $cluster2 and $size hit(s)\n";
	
	# my ($gb1, $cluster1, $query_length, $gb2,$cluster2,%hits) = @_;
	
    ##############################################################
    # Get the number of matches.
    # Calculated for each hsp as: 
    # my $num_matches = ($identity * $len ) / 100
    # and summed across all hsps.
    # This is not exact; better would be to keep track of every single base/amino acid and ask whether it was matched in some HSP.
    # But absent that, this works fairly well.
    my $num_matches = get_num_matches_bp(\%hits);

    ##############################################################
    # Total % identity across entire sequence (all HSPs).
    # Must be careful about translated queries.
    # When using %identity * length, then query length is still in bp, even for translated queries.
    # If using num matches from blast output (alignment length - num mismatches), then will be in amino acids.
    if ($query_length==0) {
    	print "Error: query length $query_length is nul for query $gb1 $cluster1 and subject $gb2 $cluster2\n";
    }
	my $identity_total = $num_matches / $query_length;

    # Print
    print OUT "$gb1 $cluster1\t";
    print OUT "$gb2 $cluster2\t";
    if($identity_total > 1) {
      print "Identity > 1 (identity = $identity_total) for $gb1 cluster$cluster1 $gb2 cluster$cluster2 \nnum matches = $num_matches, query length = $query_length\n"; 
      $identity_total = 1;
    }
    print OUT "$identity_total\t$num_matches\t$query_length\n";
}

######### END OF MODIFIED PART ###########
