#!/usr/bin/perl

# 4_generate_covs.pl
#
# Jay Moore, John Innes Centre, Norwich
# jonathan.moore@jic.ac.uk
#
# 09/10/2017
#
# Version 1
#
# Change log
#

# Summary of functions
#
# Done:
#
# Accept the project ID from command line
# Read in the project_id_sample_metadata.txt file from the 3-alignments/project_id directory
# Store the metadata in a convenience hash
# Scan the alignment folders for each sample in the project
# For each GFF file in the single-c and windows subdirectories for the sample:
# 	Parse the GFF file
# 	Replace column 3 (context ID) with an informative identifier:
# 		project_id, context, generation, line, rep 
# 	Recapitulate the GFF file structure in a new target subdirectory
#	Generate a Bismark .cov format representation of the data
#
# To do:
#


# grab the project ID from command line
my $project=$ARGV[0];

# set up a new subfolder for the re-annotateed GFF files
my $new_target = 'annotated';

# Does the alignment directory exist?
if (-e '3-alignments/'.$project) {
	if (-e '3-alignments/'.$project.'/'.$new_target) {
		print "3-alignments/$project/$new_target already exists.\n";
	} else {
		# set up convenience hash for the metadata
		my %samples;

		# parse the metadata file and grab the generation, line and rep
		open METADATA, '3-alignments/'.$project.'/'.$project.'_metadata.tsv' or die('3-alignments/'.$project.'/'.$project.'_metadata.tsv not found');
		my $line_no=0;
		while (<METADATA>) {
			my $line=$_;
			if ($line_no > 0) {
				my @sample=split("\t",$line);
				$samples{$sample[4]} = $project.'_'.$sample[1].'_'.$sample[2].'_'.$sample[3].'_'.$sample[4];
				#print $sample[4].' '.$samples{$sample[4]}."\n";
			}
			$line_no++;
		}
		close METADATA;

		# Create $new_target folder to receive outputs
		mkdir('3-alignments/'.$project.'/'.$new_target);
		
		my @subfolders;
		push @subfolders,'single-c';
		push @subfolders, 'windows';
		
		foreach my $sample(sort keys %samples) {
			mkdir('3-alignments/'.$project.'/'.$new_target.'/'.$sample.'_TAIR10');
		
			foreach my $subfolder(@subfolders) {
	
				#print $sample."\n";
				print '3-alignments/'.$project.'/'.$sample.'_TAIR10/'.$subfolder."\n";

				# Get all the GFF files for the sample for the subfolder
				opendir(DIR, '3-alignments/'.$project.'/'.$sample.'_TAIR10/'.$subfolder) or die $!;
				@files = grep { /\.gff$/ } readdir(DIR);
				closedir(DIR);
				#@list = glob("3-alignments/'.$project.'/'.$subfolder/*.GFF");
				mkdir('3-alignments/'.$project.'/'.$new_target.'/'.$sample.'_TAIR10/'.$subfolder);
				foreach my $file(@files) {
					print "$file\n";
					# Parse the GFF file and replace column 3 with a sample/context descriptive ID from the hash we prepared earlier
					open GFF, '3-alignments/'.$project.'/'.$sample.'_TAIR10/'.$subfolder.'/'.$file;
					open GFF2, '>3-alignments/'.$project.'/'.$new_target.'/'.$sample.'_TAIR10/'.$subfolder.'/'.$file;
					open COV, '>3-alignments/'.$project.'/'.$new_target.'/'.$sample.'_TAIR10/'.$subfolder.'/'.$file.'.cov';
					while (<GFF>) {
						my $line = $_;
						my @bits = split("\t",$line);
						$bits[2] = $samples{$sample}.'_'.$bits[2];
						print GFF2 join("\t",@bits);
						my @count_bits = split(';',$bits[8]);
						my @c_bits=split('=',$count_bits[0]);
						my @t_bits=split('=',$count_bits[1]);
						chomp $t_bits[1];
						if (substr($t_bits[1],-1) eq '*') {$t_bits[1]=substr($t_bits[1],0,length($t_bits[1])-1)}
						print COV $bits[0]."\t".$bits[3]."\t".$bits[4]."\t".$bits[5]."\t".$c_bits[1]."\t".$t_bits[1]."\n";
					}
					close GFF;
					close GFF2;
				}
			}		
		}
	}
} else {
	print "3-alignments/$project not found.\n";
}

