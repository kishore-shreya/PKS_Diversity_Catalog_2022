#!/usr/bin/perl -w
# Copyright (C) 2019 Maureen Hillenmeyer, Aleksandra Nivina

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

# This script was modified relative to the initial one.
# Before, there were two scripts: one that found redundancies, and another that chained them.
# Now, there's only one script that finds redundancies and chains them, in a smarter way.
# The resulting files are two lists:
# - list of non-redundant clusters, with 1 cluster/line
# - list of redundant clusters, where the first cluster/line is the "main" record, and the others are redundant
# These lists are stored in the same place as the input file with similarity scores.


use strict;

die "Usage: $0 infile cutoff\n" unless @ARGV >= 2;

my $file = shift @ARGV;
my $cutoff= shift @ARGV;

open(IN, $file);

my $outfile1 = "$file.redundant_chained_new.cutoff$cutoff.txt";
my $outfile2 = "$file.nonredundant_chained_new.cutoff$cutoff.txt";
#my $outfile3 = "$file.SequenceSimilar.cutoff$cutoff.txt"; # do this in a separate script; chaining
open(OUT1, ">$outfile1");
open(OUT2, ">$outfile2");
#open(OUT3, ">$outfile3");


# First create a hash with all similarity scores, rounded to the 3rd decimal
my %similarities;
while(<IN>) {
  chomp;
  my ($a, $b, $score, @rest) = split("\t", $_);
  next if ($a eq $b);
  $similarities{$a}{$b} = sprintf("%.3f", $score);
}


# Second, compare similarities both ways (A vs B and B vs A).
# If one of them is >cutoff: append the one that was the query to the one that was the subject.
# Ex. if A vs B is 0.99, and B vs A is 0.65, then consider B as "main" entry and A as redundant to B.
# This takes care of the situation when A is a subset of B.

my %redundant = ();
my @nonredundant = ();

# iterate over queries
foreach my $a(sort keys %similarities) {
	
	# iterate over subjects
	foreach my $b(sort keys %similarities) {
		next if ("$a" eq "$b");
		if (exists $similarities{$a}{$b}) {

		
			# check if similarity of query to subject higher than cutoff
			if ($similarities{$a}{$b}>=$cutoff) {
				
				# if so, should record this as redundant (in principle, $redundant{$b}=$a)
				
				# do nothing if this redundancy was already recorded in the other sense ($b redundant to $a already recorded, now $a redundant to $b)
					
				if (exists $redundant{$a} and check_existence_in_list($b,$redundant{$a}) == 1) {
					next;
				}
				
				# have to add redundancy:
				# either $redundant{$b}=$a if $b is non-redundant,
				# or $redundant{$d}=$a if $b is redundant to $d
				# and also add all potential records that are redundant to $a, if they have not already been saved as redundant to $b or $d
				else {				
					my $d=chained_redundancy($b,%redundant);
					# redundancies that have been recorded for $d			
					# redundancies that have been recorded for $a
					my $previous_chain= "";
					if (exists $redundant{$d}){$previous_chain= $redundant{$d};}

					my $redundant_chain= "";
					if (exists $redundant{$a}){$redundant_chain= $redundant{$a};}	
					
					# unless the query and what it's redundant to are identical, record this redundancy,
					# add all the redundancies of the query to the record that it's redundant to,
					# and delete the query as the main entry of redundancies
					unless ("$a" eq "$d") {
						# if some records redundant to $a have already been saved as redundant to $d, then only add unique ones
						# print "\n$d : $previous_chain + $a + ";
						my $additional_chain=additional_redundancy($a,$previous_chain,$redundant_chain);				
						$redundant{$d}="$previous_chain\t$additional_chain";
						delete 	$redundant{$a};
					}
				}	
			}
		}
	}
}


OUTERLOOP: foreach my $a(sort keys %similarities) {
	# if $a among keys of %redundant, then it's redundant; pass no next query
	if (exists ($redundant{$a})) {
		next OUTERLOOP;
	}
	
	# if $a is among values of %redundant, then it's redundant; pass no next query
	foreach my $b(sort keys %redundant) {
		if (check_existence_in_list($a,$redundant{$b}) == 1) {
			next OUTERLOOP;
		}

	}
	
	# otherwise, record $a as non-redundant
	push (@nonredundant, $a);
}

# Recording redundant into the OUT1 file
for my $main_redund (sort keys %redundant) {	
	my @other_redund=split ("\t",$redundant{$main_redund});	
    print OUT1 "$main_redund\t";
    print OUT1 join ("\t", @other_redund);
    print OUT1 "\n";
}

# Recording non-redundant into the OUT1 file
for my $nonredund (@nonredundant) {
	print OUT2 "$nonredund\n";
}



# my $nr =length (@nonredundant);
my $redund= keys %redundant;
print "\n";
print "non-redundant = $#nonredundant\n";
print "redundant sets = $redund\n";

#print "\noutput to \n$outfile1\nand\n$outfile2\nand\n$outfile3\n\n";
print "\noutput to \n$outfile1\nand\n$outfile2\n\n\n";

exit;


close(IN);
close(OUT1);
close(OUT2);

exit;



sub chained_redundancy {
	my ($b, %hash) = @_;
	
	# iterating over $d: clusters that are already recorded as having redundancies
	foreach my $d (sort keys %hash) {
		
		# if $b is recorded as redundant to some other $d, then we return this $d, because $a should be redundant to it
		if (check_existence_in_list($b,$hash{$d}) == 1) {
			return $d;
		}
	}
	
	# otherwise, $b isn't redundant to anything, then we return this $b, because $a should be redundant to it
	return $b;
}

sub additional_redundancy {
	my $a=shift;
	my $previous_chain=shift;
	my $redundant_chain=shift;
	my @previous_records=split("\t",$previous_chain);
	my @additional_records=split("\t",$redundant_chain);
	my @unique_additional_records=();

	# if $a not is already present in the previous records, add it
	my $a_presence=0;
	foreach my $previous_record (@previous_records) {
		if ($a eq $previous_record) {
			$a_presence=1;
			last;
		}	
	} 
	if (not $a_presence) {
		push (@unique_additional_records,$a);
	}
	
	# if a record not is already present in the previous records, add it
	OUTLOOP: foreach my $record (@additional_records) {
		INLOOP: foreach my $previous_record (@previous_records) {
			if ($record eq $previous_record) {
				next OUTLOOP;
			}
		}
		push (@unique_additional_records,$record);
	}
	my $additional_chain=join ("\t", @unique_additional_records);
	# print "$additional_chain\n";
	return $additional_chain;
	
	
}

sub check_existence_in_list {
    my $check_existence_accession = shift;
    my $string_lst_accession = shift;
    my @lst_accession=split("\t", $string_lst_accession);
    my $result = 0;
    foreach my $accession (@lst_accession){
        if($accession eq $check_existence_accession){
            $result++;
        }
    }
    return $result;
}
