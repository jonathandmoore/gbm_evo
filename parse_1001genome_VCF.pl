#!/usr/bin/perl

# parse_1001genome_VCF.pl
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 08/11/2019
#
# Version 1
#
# Change log
#

# Summary of functions
#
# Done:
#
# Read .bed file of CG sites
# Stream .VCF file on STDIN
# Parse the file identifying cytosines and CG dinucleotides
# Output separate lists, one for each accession, of CG sites with mutations and of cytosines with mutations
#
# To do:
#

my $cg_sites_file = '/jic/scratch/groups/Daniel-Zilberman/Reference_data/Genomes/Arabidopsis/TAIR10/TAIR10_Chr.all.CG_sites.tsv';
# make a hash to put the CG sites in
my %cg_sites;
open(my $cg_sites_fh, '<:encoding(UTF-8)', $cg_sites_file) or die "Could not open file '$cg_sites_file' $!";
while (<$cg_sites_fh>) {
  my $line = $_;
  chomp $line;
  my @this_site = split "\t", $line;
  # create an entry in the hash for the CG site pasting together chrom and locus1 as the key
  $cg_sites{$this_site[0].$this_site[1]}=1;
}
close $cg_sites_fh;

#foreach my $cg_site(sort keys %cg_sites) {
#  print "$cg_site\n";
#}

my %samples;
while (<>) {
  my $line = $_;
  chomp $line;
  if (substr($line,0,1) eq '#') {
    if (substr($line,0,2) eq '##') {
      # skip comments
	} else {
	  # it's the header line with all the accession IDs - grab them for later
	  my @header = split '\t', $line;
	  for (my $sample=9; $sample<=$#header; $sample++) {
	    $samples{$sample} = $header[$sample];
	  }
	  foreach my $sample(sort keys %samples) {
	  #  print $sample." ".$samples{$sample}."\n";
	    system('rm '.$samples{$sample}.'_CG_vars.txt');
	    system('rm '.$samples{$sample}.'_C_vars.txt');
	  }
	}
  } else {
    # we have finished headers and are now doing content
	my @entries = split '\t', $line;
    # make list of lines mutated at this site - start with an empty array
	my @mutated_lines = ();
	print $entries[0].' '.$entries[1].' '.$entries[2].' '.$entries[3].' '.$entries[4].' '.$entries[5].' '.$entries[6].' '.$entries[7].' '.$entries[8]."\n";

    for (my $sample=9; $sample<=$#entries; $sample++) {
	  # parse the entry of a single accession
	  my @entry = split ':', $entries[$sample];
      #print "entry[0] ".$entry[0]."\n";
	  if (($entry[0] ne '0|0') & ($entry[0] ne './.')) {
	    # if the entry is neither ref|ref nor missing data, it is a mutation - log the relevant sample column number to the array
        push @mutated_lines, $sample;
		print $samples{$sample}." ".$entry[0]."\n";
	  }
	}

    # identify the sequence identified as the reference version
	my $refseq = $entries[3];

	# go through the mutated lines writing out the site info to the relevant files
    foreach my $mutant (@mutated_lines) {
	  #print "hello $mutant\n";
	  
	  #if ((length $samples{$mutant}) == 0) {
	    #print "hello $mutant\n";
	  #}

	  # do any of the sites in the reference version overlap a CG site?
	  my $overlaps_cg = 0;
	  #print length $refseq;
      for (my $this_site=0; $this_site<(length $refseq); $this_site++) {
	    # if either the site or its neighbour are the cytosine of a CG site
        if (defined $cg_sites{'Chr'.$entries[0].$entries[1]}) {
		  if ($overlaps_cg != $entries[1]) {
            $overlaps_cg = $entries[1];
    		open(my $fh, '>>', $samples{$mutant}.'_CG_vars.txt') or die "Could not write to $samples{$mutant}_CG_vars.txt $!";
	    	say $fh 'Chr'.$entries[0]."\t".$entries[1];
		    close $fh;
		  }
        } elsif (defined $cg_sites{'Chr'.$entries[0].($entries[1]-1)}) {
          # has this CG site already been reported?
		  # There is still a chance that the same CG site gets reported twice, if both nucleotides in the site have separate mutations. In that case the site will get reported twice in any genotypes which have a non-reference allele at both sites
		  if ($overlaps_cg != ($entries[1]-1)) {
            $overlaps_cg = $entries[1]-1;
    		open(my $fh, '>>', $samples{$mutant}.'_CG_vars.txt') or die "Could not write to $samples{$mutant}_CG_vars.txt $!";
	    	say $fh 'Chr'.$entries[0]."\t".($entries[1]-1);
		    close $fh;
		  }
		}
      }
	  #if ($overlaps_cg > 0) {
        # output the relevant CG site locus to the files for the samples which vary from the reference
		#print 'hello';
		#open(my $fh, '>>', $samples{$mutant}.'_CG_vars.txt') or die "Could not write to $samples{$mutant}_CG_vars.txt $!";
		#say $fh 'Chr'.$entries[0]."\t".$entries[1];
		#close $fh;
      #}
      
	  # are any of the sites in the reference version a C (or a G) in the reference?
	  my $overlaps_c = 0;
      for (my $this_site=0; $this_site<(length $refseq); $this_site++) {
        if ((substr($refseq,$this_site,1) eq 'C') | (substr($refseq,$this_site,1) eq 'G')) {
		  $overlaps_c = 1;
		  $this_subsite = $entries[1]+$this_site;
          # output the relevant C/G site locus to the files for the samples which vary from the reference
	  	  open(my $fh, '>>', $samples{$mutant}.'_C_vars.txt') or die "Could not write to $samples{$mutant}_C_vars.txt $!";
		  say $fh 'Chr'.$entries[0]."\t".$this_subsite;
		  close $fh;
	    }
	  }
	  #if ($overlaps_c == 1) {
        # output the relevant C site locus to the files for the samples which vary from the reference
		#open(my $fh, '>>', $samples{$mutant}.'_C_vars.txt') or die "Could not write to $samples{$mutant}_C_vars.txt $!";
		#say $fh 'Chr'.$entries[0]."\t".$entries[1];
		#close $fh;
      #}
	}
  }
}

