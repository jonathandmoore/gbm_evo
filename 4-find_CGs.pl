#!/usr/bin/perl

# 4-find_CGs.pl
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 18/10/2017
#
# Version 1
#
# Change log
#

# Summary of functions
#
# Done:
#
# Accept a FASTA filepath from command line
# Parse the file identifying CG dinucleotides
# Output a list of tab-separated Chromosome IDs and start/end loci of CG dinucelotides
#
# To do:
#


# grab the FASTA filename from command line
my $fasta=$ARGV[0];
# This is the pattern to search for
my $pattern = 'CG';
# This is the adjustment to base position (0 or 1 base)
my $adjustment = 1;

# Does the alignment directory exist?
if (-e $fasta) {
	# load the FASTA file into a hash, one entry per sequence, with no return characters
	my %fasta;
	open FASTA, $fasta;
	my $current_sequence;
	while (<FASTA>) {
		my $line=$_;
		chomp $line;
		if (substr($line,0,1) eq '>') {
			$current_sequence = substr $line, 1;
			my @current_sequence = split ' ',$current_sequence;
			$fasta{$current_sequence[0]} = '';
			$current_sequence = $current_sequence[0];
			#print "Loading $current_sequence\n";
		} else {
			$fasta{$current_sequence} .= $line;
		}
	}
	close FASTA;
	
	# Go through the entries in the hash outputting CG site loci
	foreach my $seq(sort keys %fasta) {
	#print $seq.' '.length($fasta{$seq})."\n";
	for (my $locus=0; $locus<length($fasta{$seq}); $locus++) {
			if (uc(substr($fasta{$seq},$locus,2)) eq $pattern) {
				my $start = $locus+$adjustment;
				my $end = $locus+$adjustment+1;
				print $seq."\t".$start."\t".$end."\n";
			}
		}
	}
} else {
	print "$fasta not found.\n";
}

