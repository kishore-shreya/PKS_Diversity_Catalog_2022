#!/usr/bin/perl -w
# NOTE: THIS contains an error in the overlap function
# See ~/TypeIICatalog/findClusteredKSCLF.pl for corrected overlap function
use strict;
use Bio::SearchIO;
use List::Util qw(min max);

die "Usage: $0 blast_filename\n" unless @ARGV;

my $infile = shift @ARGV;
my $distCutoff = 20000;# for clustering nearby KSs
my $tooCloseCutoff = 3000; # for filtering KSs that are too close together, not separate modules

my %clusterMin;
my %clusterMax;
my %genbank_hspNum_clusterNum; # Reassigned cluster num. key = genbank, key = hspNum, val = clusterNum
my %genbank_clusterNum_hspNums; # key = genbank, key = clusterNum, keys = hspNums
my %genbank_desc = ();
my %genbank_length = ();
my %tooClose = (); #{$genbank}{$hspX}{$hspY}

my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $infile);


warn "Iterating blast results...\n";

while( my $result = $in->next_result ) {

  while( my $hit = $result->next_hit ) {

    my $hitname = $hit->name;
    my $hitdesc = $hit->description;
    warn "$hitdesc";
    #next unless $hitdesc =~ /Toxo/; ####################### debug

    my $hitlen = $hit->length;
    $genbank_desc{$hitname} = $hitdesc;
    $genbank_length{$hitname} = $hitlen;
    my @hsps = $hit->hsps;

    my $hspXNum = 0;
    foreach my $hsp(@hsps) {
      my $frame = $hsp->hit->frame;
      my $strand= $hsp->strand('hit');
      my $start = $hsp->start('hit');
      my $end = $hsp->end('hit');
      $genbank_clusterNum_hspNums{$hitname}{$hspXNum}{$hspXNum} = 1; # 
      $genbank_hspNum_clusterNum{$hitname}{$hspXNum} = $hspXNum; # at first, clusterID is hspID
      $clusterMin{$hitname}{$hspXNum} = $start; #unless strand = -1? no ok
      $clusterMax{$hitname}{$hspXNum} = $end; #unless strand = -1? no ok
      my $hspYNum = 0;

      # Iterate HSPs again, calculate distance between HSP pairs
      foreach my $hspY(@hsps) {
        if ($hsp == $hspY) {$hspYNum++; next;}
 	my $strandY= $hspY->strand('hit');
        my $startY = $hspY->start('hit');
        my $endY = $hspY->end('hit');
	my $minX = min($start, $end);
	my $minY = min($startY, $endY);
	my $maxX = max($start, $end);
	my $maxY = max($startY, $endY);
	my $dist = abs($minX - $minY);
	my $cluster = "nocluster";

	# debug
	#next unless ($hitname =~ /AAQM02000016/);
	#warn "$hitname\t$hspXNum\t$hspYNum\t$dist\n";

	# If nearby, then merge
	if($dist < $distCutoff) {
  	  if($dist < $tooCloseCutoff) {
		#need to mark them as one so they don't get counted as separate
		$tooClose{$hitname}{$hspXNum}{$hspYNum} = 1;
		#warn "too close $hitname $hspXNum $hspYNum\n";
	  }

          $genbank_clusterNum_hspNums{$hitname}{$hspXNum}{$hspYNum} = 1; #

	  # Assign this hspY to the earliest hspX that it can find -- trace back the chain
          my $hspsX = $genbank_clusterNum_hspNums{$hitname}{$hspXNum};
	  my @hspsX = sort keys %$hspsX;
	  my $earliestHspXNum = min(@hspsX);
          $genbank_hspNum_clusterNum{$hitname}{$hspYNum} = $earliestHspXNum; # key = genbank, key = hspNum, val = clusterNum

	  # Reset cluster min/max if necessary
	  my $currentMin = min($minX, $minY);
	  my $currentMax = max($maxX, $maxY);
	  my $clusterMin = $clusterMin{$hitname}{$earliestHspXNum};
	  my $clusterMax = $clusterMax{$hitname}{$earliestHspXNum};
	  if($currentMin < $clusterMin) {
		$clusterMin{$hitname}{$earliestHspXNum} = $currentMin;
	  }
	  if($currentMax > $clusterMax) {
		$clusterMax{$hitname}{$earliestHspXNum} = $currentMax;
	  }

	  $cluster = "cluster";
	}
	$hspYNum++;
      }
    $hspXNum++;
    }
  }
}



# Iterate clusters, merging overlapping ones

warn "Iterating clusters, merging overlapping ones...\n";

foreach my $genbank (sort keys %genbank_clusterNum_hspNums) {
  my $clusterNums = $genbank_clusterNum_hspNums{$genbank};

  # Iterate HSPs
  foreach my $clusterNumX(sort {$a <=> $b} keys %$clusterNums) {
    my $hspsX = $genbank_clusterNum_hspNums{$genbank}{$clusterNumX};
    my @hspsX = sort {$a <=> $b} keys %$hspsX;
    my $numHSPsX = @hspsX;
    my $clusterMinX = $clusterMin{$genbank}{$clusterNumX};
    my $clusterMaxX = $clusterMax{$genbank}{$clusterNumX};

    #print "$genbank\t$clusterNumX\t@hspsX\t$clusterMinX\t$clusterMaxX\n";

    my $earliestOverlapNumSeen = 1000000000;

    # Iterate HSPs again, looking for overlaps to merge
    foreach my $clusterNumY (sort {$a <=> $b} keys %$clusterNums) {
      next if $clusterNumX == $clusterNumY;
      my $hspsY = $genbank_clusterNum_hspNums{$genbank}{$clusterNumY};
      my @hspsY = sort {$a <=> $b} keys %$hspsY;
      my $numHSPsY = @hspsY;
      my $clusterMinY = $clusterMin{$genbank}{$clusterNumY};
      my $clusterMaxY = $clusterMax{$genbank}{$clusterNumY};

      if(overlap($clusterMinX, $clusterMaxX, $clusterMinY, $clusterMaxY)) {
        #print "  OVERLAP : $clusterNumY\t@hspsY\t$clusterMinY\t$clusterMaxY\n";

	# Reassign
	$genbank_clusterNum_hspNums{$genbank}{$clusterNumX}{$clusterNumY} = 1; #
        $earliestOverlapNumSeen = min($earliestOverlapNumSeen, @hspsX, @hspsY);
	my $earliestHsp = min(@hspsX, @hspsY, $earliestOverlapNumSeen);
	#print "  Earliest HSP = $earliestHsp\n";
	$genbank_clusterNum_hspNums{$genbank}{$earliestHsp}{$clusterNumY} = 1; #
	$genbank_clusterNum_hspNums{$genbank}{$earliestHsp}{$clusterNumX} = 1; # necessary?  already there?
        $genbank_hspNum_clusterNum{$genbank}{$clusterNumY} = $earliestHsp; # key = genbank, key = hspNum, val = clusterNum
        $genbank_hspNum_clusterNum{$genbank}{$clusterNumX} = $earliestHsp; # key = genbank, key = hspNum, val = clusterNum

	# Reset cluster min/max if necessary
	if($clusterMinY < $clusterMin{$genbank}{$clusterNumX}) { $clusterMin{$genbank}{$clusterNumX} = $clusterMinY; }
	if($clusterMaxY > $clusterMax{$genbank}{$clusterNumX}) { $clusterMax{$genbank}{$clusterNumX} = $clusterMaxY; }
	if($clusterMinY < $clusterMin{$genbank}{$earliestHsp}) {$clusterMin{$genbank}{$earliestHsp} = $clusterMinY;}
	if($clusterMaxY > $clusterMax{$genbank}{$earliestHsp}) {$clusterMax{$genbank}{$earliestHsp} = $clusterMaxY;}
	my $hspsNew = $genbank_clusterNum_hspNums{$genbank}{$earliestHsp};
    	my @hspsNew = sort {$a <=> $b} keys %$hspsNew;
      }
    }
    my $reassigned = $genbank_hspNum_clusterNum{$genbank}{$clusterNumX};
  }  
}

# Gather unique cluster IDs
my %uniqueClusters; #key = genbank, key = unique clusterid

foreach my $genbank (sort keys %genbank_clusterNum_hspNums) {
  my $clusterNums = $genbank_clusterNum_hspNums{$genbank};
  foreach my $clusterNumX(sort {$a <=> $b} keys %$clusterNums) {
    my $hspsX = $genbank_clusterNum_hspNums{$genbank}{$clusterNumX};
    my @hspsX = sort {$a <=> $b} keys %$hspsX;
    my $numHSPsX = @hspsX;
    my $clusterMinX = $clusterMin{$genbank}{$clusterNumX};
    my $clusterMaxX = $clusterMax{$genbank}{$clusterNumX};
  
    my $reassignedCluster = $genbank_hspNum_clusterNum{$genbank}{$clusterNumX};
    $uniqueClusters{$genbank}{$reassignedCluster} = 1; 

    #print "$genbank\t$clusterNumX\treassigned: $reassignedCluster\t@hspsX\t$clusterMinX\t$clusterMaxX\n";

  }
}


foreach my $genbank (sort keys %uniqueClusters) {
  my $desc = $genbank_desc{$genbank};
  my $length = $genbank_length{$genbank};
  my $clusters = $uniqueClusters{$genbank};
  foreach my $cluster(sort {$a <=> $b} keys %$clusters) {
    my $hsps = $genbank_clusterNum_hspNums{$genbank}{$cluster};
    my @hsps = sort {$a <=> $b} keys %$hsps;

    #Collapse num of HSPs by checking which ones were too close (<3kb)
    my @hspsFiltered = filterHspsTooClose($genbank, \@hsps);

    my $numHsps = @hspsFiltered;
    my $clusterMinX = $clusterMin{$genbank}{$cluster};
    my $clusterMaxX = $clusterMax{$genbank}{$cluster};
    #print "$genbank\t$desc\t$length\t$cluster\t$numHsps\t@hsps\t$clusterMinX\t$clusterMaxX\n";
    print "$genbank\t$desc\t$length\t$cluster\t$numHsps\t@hsps\t@hspsFiltered\t$clusterMinX\t$clusterMaxX\n";
  }
}




exit;



###########################################
# Subroutines

sub filterHspsTooClose{
  my $genbank = shift;
  my($hsps) = shift;
  my @hsps = @$hsps; # already sorted
  my @hspsFiltered = ();
  my %seen = ();
  #my %toss = ();
  my %combine = ();
  foreach my $hspX (@hsps) {
    foreach my $hspY(@hsps) {
      my $tooClose = $tooClose{$genbank}{$hspX}{$hspY};
      if($tooClose) {
	#$toss{$hspY} = 1;
	#warn "Tossing $hspY\n";
	$combine{$hspX}{$hspY} = 1;
      }
    }
  }

  foreach my $hspX(@hsps) {
    #unless ($toss{$hspX}) { push(@hspsFiltered, $hspX); }
    next if $seen{$hspX};
    if($combine{$hspX}) {
      my $combine_str = "$hspX";
      # Get the hsps to combine X with
      my $hspsY = $combine{$hspX};
      foreach my $hspY(keys %$hspsY) {
	$combine_str .= ",$hspY";
        $seen{$hspY} = 1;
      }
      push(@hspsFiltered, $combine_str);
    }
    else {
      push(@hspsFiltered, $hspX);
    }
  }
  return (@hspsFiltered);
}


sub overlap{ ######### note this has a bug.  Doesn't allow for the case in which Y is entirely contained within X (X max and min are outside of Y's max and min)
  my ($clusterMinX, $clusterMaxX, $clusterMinY, $clusterMaxY) = @_;
  #print "$clusterMinX, $clusterMaxX, $clusterMinY, $clusterMaxY\n";
  if ($clusterMaxX >= $clusterMinY && $clusterMaxX <= $clusterMaxY ||
      $clusterMinX >= $clusterMinY && $clusterMinX <= $clusterMaxY) {
	return 1;
  }
  return 0;
}


