#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(min max sum);
use FindBin;
use lib "$FindBin::Bin/lib";
use feature 'say';
#use DZUtil qw/read_conf/;


my $DATA_HANDLE    = 'ARGV';
my $gff_annotation = q{};
my $bin_width      = 100;
my $distance       = 5000;
my $stop_flag      = 2;
my $stop_distance  = 1500;
my $three_prime    = 0;
my $five_prime     = 0;
my $attribute_id   = 'ID';
my $output         = '-';
my $debug;
my $zero_flag_region;

# gff-annotation = /path/to/annotation.gff
# bin-width      = 100
# distance       = 5000
# stop-flag      = 6 
# stop-distance  = 1500
# end            = 5
# extract-id     = ID

#my %conf = read_conf();
#$gff_annotation = $conf{'gff-annotation'} // $gff_annotation;
#$bin_width     = $conf{'bin-width'}     // $bin_width;
#$distance      = $conf{'distance'}      // $distance;
#$stop_flag     = $conf{'stop-flag'}     // $stop_flag;;
#$stop_distance = $conf{'stop-distance'} // $stop_distance;
#$attribute_id  = $conf{'extract-id'}    // $attribute_id;
#if (exists $conf{end}){
#    if (3 == $conf{end}){
#        $three_prime = 1;
#    } elsif (5 == $conf{end}){
#        $five_prime = 1;
#    }
#}

# Grabs and parses command line options
my $result = GetOptions(
    'gff-annotation|g=s' => \$gff_annotation,
    'bin-width|b=i'      => \$bin_width,
    'distance|d=i'       => \$distance,
    'stop-flag|s=i'      => \$stop_flag,
    'stop-distance|k=i'  => \$stop_distance,
    'three-prime|3'      => \$three_prime,
    'five-prime|5'       => \$five_prime,
    'extract-id|x=s'     => \$attribute_id,
    'output|o=s'         => \$output,
    'debug'              => \$debug,
    'zero'               => \$zero_flag_region,
    'no-skip'            => \(my $noskip),
    'verbose|v'          => sub { use diagnostics; },
    'quiet|q'            => sub { no warnings; },
    'help|h'             => sub { pod2usage( -verbose => 1 ); },
    'manual|m'           => sub { pod2usage( -verbose => 2 ); }
);

say STDERR "gff_annotation $gff_annotation";
say STDERR "bin_width      $bin_width     ";
say STDERR "distance       $distance      ";
say STDERR "stop_flag      $stop_flag     ";
say STDERR "stop_distance  $stop_distance ";
say STDERR "three_prime    $three_prime   ";
say STDERR "five_prime     $five_prime    ";
say STDERR "attribute_id   $attribute_id  ";
say STDERR "output         $output        ";

# Add additional stop strategies here and implement the corresponding function
my $stop_flag_dispatch = {
    0 => \&stop_flag_0,    # implemented; passed tests
    1 => \&stop_flag_1,
    2 => \&stop_flag_2,    # implemented; check embedded loci
    3 => \&stop_flag_3,
    4 => \&stop_flag_4,
    5 => \&stop_flag_5,
    6 => \&stop_flag_6,    # implemented: see #2
};

# Check required command line parameters
pod2usage( -verbose => 1 )
    unless (@ARGV
    and $result
    and $gff_annotation
    and ( $three_prime xor $five_prime )
    and ( $stop_flag != 6 or ( $stop_flag == 6 and $stop_distance ) )
    and exists $stop_flag_dispatch->{$stop_flag},
    and $bin_width > 0);

if ($output ne '-') {
    open my $USER_OUT, '>', $output or croak "Can't read $output: $!";
    select $USER_OUT;
}

my %flag_parameters = (
    -three_prime   => $three_prime,
    -five_prime    => $five_prime,
    -distance      => $distance,
    -stop_flag     => $stop_flag,
    -stop_distance => $stop_distance,
);

# index and offset the loci annotations data based on stop flag
# { seqname => { 
#      start_coord => 
#           [flag_start, flag_end, strand, attributes, 
#            original_start, original_end]}}
# attribute = locus id
my $index_annotation = index_gff_annotation( $gff_annotation, $attribute_id );

my $annotation = offset_gff_annotation(
    $index_annotation,
    $stop_flag_dispatch->{$stop_flag},
    \%flag_parameters
);

#say STDERR Dumper $annotation;
say STDERR scalar map { values %$_ } values %$index_annotation;

#######################################################################
# locus->coord-matrix lookup 

my %gene_info;
my %gene_seen;

while (my ($seqname,$start2coords) = each %$index_annotation) {
    while (my ($start,$coords) = each %$start2coords) {
        $gene_seen{$coords->[3]} = 0;
    }
}
while (my ($seqname,$start2coords) = each %$annotation) {
    while (my ($start,$coords) = each %$start2coords) {
        $gene_info{$coords->[3]} = $coords;
    }
}



#say STDERR Dumper $annotation;

my $num_bins           = 2 * int( $distance / $bin_width );
my %genes              = (); # hash of loci id hashes with scores/strand keys
my %sorted_annotations = (); # cache for sorted GFF annotations
my $gff_iterator       = make_gff_iterator( $ARGV[0], \&gff_read );

COORD:
while ( my $gff_line = $gff_iterator->() ) {

    # Step 1: Parsing GFF line
    # &gff_read returns [] for GFF comments, invalid lines, etc.
    next COORD unless ref $gff_line eq 'HASH';

    next COORD if $gff_line->{score} eq q{.};

    # Step 2: Memoization for sorting annoation
    # has this been done before? use it : otherwise cache it
    unless ( exists $sorted_annotations{ $gff_line->{seqname} } ) {

        my @loci_coords   = keys %{ $annotation->{ $gff_line->{seqname} } };
        my @sorted_coords = sort { $a <=> $b } @loci_coords;

        # map to values [start, end, strand, locus id]
        $sorted_annotations{ $gff_line->{seqname} } = [
            map { $annotation->{ $gff_line->{seqname} }{$_} }
                @sorted_coords
        ];
    }

    # Step 3: Searching loci where this score fits into
    # the look-up key is a range reference, brs searches an array of range references
    my $locus = binary_range_search(
        [ $gff_line->{start}, $gff_line->{end} ],
        $sorted_annotations{ $gff_line->{seqname} },
    ) || next COORD;

    # Step 4: Reverse search if needed
    # orientation of search: 5'->3' or 5'<-3' (reverse)
    my $reverse
        = ( $five_prime and $locus->[2] eq q{-} )
            || ( $three_prime and $locus->[2] eq q{+} )
                || 0;

    # Step 5: Locate bin of width bin_width where this score fits into
    # the bin index is the locus 5' or 3' coordinate
    # minus the current start coordinate
    # minus bin_width * (num_bins / 2)
    # divided by the negative bin width.
    # note: $reverse + 4 accesses original start or end coordinate for this locus
    my $index = int(
        ($locus->[ $reverse + 4 ]          # original start
         - $bin_width * ( $num_bins / 2 )  # new start. is just $distance, no?
         - $gff_line->{start}              # result diff between $gff_line location and new start. negative
        ) / -$bin_width
    );
    #print("hello $locus->[ $reverse + 4 ] $bin_width * ( $num_bins / 2 ) $gff_line->{start}  hello\n");

    # Step 6: Weight score by amount of overlap
    my $score = weight_score( $locus, $gff_line );

    # Step 7: Save score into its appropriate bin, record its strand
    #print("hello $index hello\n");
	# Added this test because it was causing an error if $index was negative in the first bin where it saw a gene
    if (($index>=0) | ($gene_seen{$locus->[3]} == 1)) {
      push @{ $genes{ $locus->[3] }->{scores}[$index] }, $score;
  	  $gene_seen{$locus->[3]} = 1;
    }
}

#say STDERR Dumper \%genes;

# Step 8: Iterate over scores for each gene and print those
#say Dumper \%genes;
my $scores_iterator = make_scores_iterator( \%genes );
while ( my ( $locus, $scores_ref ) = $scores_iterator->($num_bins) ) {
    if ($debug){
        my $counter = -$distance;
        print "$locus\n";
        my @scores_copy = @$scores_ref;
        while (scalar @scores_copy){
            print "$counter\t";
            for (1 .. 10){
                my $score = shift @scores_copy;
                print "\t:$score";
                $counter+=$bin_width;
            }
            print "\n";
        }
        print "\n";
    }
    else{
        print join( "\t", $locus, @{$scores_ref} ), "\n";
    }
}

## done

#==================================================================
# only assign partial scores according to amount of overlap

sub weight_score {
    my ( $locus, $gff_line ) = @_;

    my $bin_low  = max( $locus->[0], $gff_line->{start} );
    my $bin_high = min( $locus->[1], $gff_line->{end} );
    my $overlap  = abs($bin_high - $bin_low) + 1;
    my $weight   = $overlap / ( $gff_line->{end} - $gff_line->{start} + 1 );

    return $gff_line->{score} * $weight;
}

sub make_scores_iterator {
    my ($genes_ref) = @_;
    #die Dumper $genes_ref;
    my %local_genes = %$genes_ref;

    my @genes_names = $noskip ? sort keys %gene_seen : sort keys %local_genes;

    return sub {
        my ($num_bins) = @_;

        my $gene = shift @genes_names || return;
        my @scores = ();

        # if gene was not seen AND there's no entry in gene_info,
        # it means it didn't get past offset_annotation b/c it didn't lead to 
        # proper flag start/end.  Else, it can go forward even if it wasn't seen.
        if ($noskip && ! $gene_seen{$gene} && ! exists $gene_info{$gene}){
            return ($gene, [('na') x $num_bins]);
        }

        #say Dumper $gene_info{$gene};
        my $gene_start = $gene_info{$gene}[4];
        my $gene_end   = $gene_info{$gene}[5];
        my $flag_start = $gene_info{$gene}[0];
        my $flag_end   = $gene_info{$gene}[1];
        my $strand     = $gene_info{$gene}[2];
        if ($debug){
            say "gene: $gene";
            say "gene_start: $gene_start";
            say "gene_end: $gene_end";
            say "flag_start: $flag_start";
            say "flag_end: $flag_end";
            say "strand: $strand";
        }

        for my $index ( 0 .. $num_bins - 1 ) {
            # if '-', go from $num_bins -> 0
            $index =  $strand eq q{-} ? $num_bins - $index - 1 : $index;

            my $gene_forward_relcoord = 
            ( (($strand eq '+' && $five_prime) || ($strand eq '-' && $three_prime)) ?  $gene_start : $gene_end) 
            - $distance + $index * $bin_width;
            my $in_flag = $flag_start <= $gene_forward_relcoord && $gene_forward_relcoord <= $flag_end;

            say "$index, $gene_forward_relcoord, $in_flag"
                if ($debug);

            if (defined $local_genes{$gene}->{scores}[$index]){
                push @scores,
                sprintf( "%g",
                    ( sum @{ $local_genes{$gene}{scores}[$index] } )
                        / @{ $local_genes{$gene}{scores}[$index] } ) ;
            } 
            elsif($zero_flag_region && $in_flag){
                push @scores, 0;
            } 
            else{
                push @scores, 'na';
            }
        }
        return ( $gene, \@scores );
    }
}

sub binary_range_search {
    my ( $range, $ranges ) = @_;
    my ( $low, $high ) = ( 0, @{$ranges} - 1 );
    while ( $low <= $high ) {

        my $try = int( ( $low + $high ) / 2 );

        $low = $try + 1, next if $ranges->[$try][1] < $range->[0];
        $high = $try - 1, next if $ranges->[$try][0] > $range->[1];

        # range_assert ($range->[0], @{$ranges->[$try]});
        return $ranges->[$try];
    }
    return;
}

sub range_assert {
    my ( $coord, $low, $high ) = @_;

    if ( $coord < $low or $coord > $high ) {
        croak "NOT OK: ", $coord, "\t$low-$high\n";
    }
}

sub make_gff_iterator {
    my ( $gff_file_name, $gff_parser ) = @_;

    open my $GFF_FILE, '<', $gff_file_name
        or croak "Can't read $gff_file_name: $!";

    return sub { $gff_parser->( scalar <$GFF_FILE> ) };
}

sub gff_read {
    return [] if $_[0] =~ m/^
                            \s*
                            \#+
                           /mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, shift || return;

    $attribute =~ s/[\r\n]//mxg;

    return {
        'seqname'   => lc $seqname,
        'source'    => $source,
        'feature'   => $feature,
        'start'     => $start,
        'end'       => $end,
        'score'     => $score,
        'strand'    => $strand,
        'frame'     => $frame,
        'attribute' => $attribute
    };
}

#==================================================================
# for a given gff file and an attribute tag, return 
# { seqname => { start_coord => [start, end, strand, attribute}}

sub index_gff_annotation {
    my ( $gff_file_name, $attribute_id ) = @_;

    my $gff_annotation = {};
    my $gff_iterator = make_gff_iterator( $gff_file_name, \&gff_read );

 GFF:
    while ( my $gff_line = $gff_iterator->() ) {

        next GFF unless ref $gff_line eq 'HASH';

	$gff_line->{attribute} =~ s/.*
                                    $attribute_id
                                    =
                                    (\w+)
                                    .*
                                   /$1/mx
                                       if $attribute_id;

        $gff_annotation->{ $gff_line->{seqname} }{ $gff_line->{start} } = [
            $gff_line->{start},  $gff_line->{end},
            $gff_line->{strand}, $gff_line->{attribute}
        ];
    }
    return $gff_annotation;
}

#==================================================================
# like index_gff_annotation EXCEPT replace the start/end with flag start/end.
# keep original start/end in array index 4/5.
# return { seqname => { 
#             start_coord => 
#                  [flag_start, flag_end, strand, attributes, 
#                   original_start, original_end]}}
sub offset_gff_annotation {
    my ( $gff_annotation, $flag_parser, $parameters ) = @_;
    # gff_annotion: { seqname => { start_coord => [start, end, strand, attribute}}
    # flag_parser: subref which returns flag start/end
    # parameters: -three_prime -five_prime -distance -stop_flag -stop_distance (-stop_flag redundant)

    my $offset_gff_annotation = {};

    return unless ($gff_annotation and $flag_parser);

    for my $seqid ( sort keys %{$gff_annotation} ) {

        # every count of 3 genes, we offset the middle one
        # this ensures that we don't miss the first gene
        my @memory = ( [ 0, 0, q{+}, q{} ] );

        for my $start (
            sort { $a <=> $b }
                keys %{ $gff_annotation->{$seqid} }
            ) {

            push @memory, $gff_annotation->{$seqid}{$start};

            # buffer is full, offset middle gene and delete first gene
            if ( @memory == 3 ) {
                check_neighbourhood( $offset_gff_annotation, $flag_parser,
                    $parameters, $seqid, @memory );
                shift @memory;
            }
        }

        # this ensures that we don't miss the last gene
        push @memory, [ 0, 0, q{+}, q{} ];
        check_neighbourhood( $offset_gff_annotation, $flag_parser,
            $parameters, $seqid, @memory );

    }
    return $offset_gff_annotation;
}


sub check_neighbourhood {
    my ( $offset_gff_annotation, $flag_parser, $parameters, $seqid, @memory )
        = @_;

    # offset gene coordinates using appropriate flag method
    my ( $flag_start, $flag_end ) = $flag_parser->( $parameters, @memory );

    $offset_gff_annotation->{$seqid}{$flag_start} = [
        $flag_start,     $flag_end,       $memory[1]->[2],
        $memory[1]->[3], $memory[1]->[0], $memory[1]->[1]
    ]
    if $flag_start <= $flag_end;

    return 1;
}



sub stop_flag_0 {
    my ( $parameters, $previous, $current, $next ) = @_;

    return ( min_max_distance( $current, $parameters ) );
}

sub stop_flag_1 {
    croak "Not yet implemented";
}

sub stop_flag_2 {
    my ( $parameters, $previous, $current, $next ) = @_;

    my ( $minimum, $maximum ) = min_max_distance( $current, $parameters );

    # basically: go down or upstream at most -distance, but clip so that you 
    # don't go beyond:
    # 1) the other end of itself,
    # 2) into another gene
    if (   $parameters->{-five_prime} and $current->[2] eq q{-}  
        or $parameters->{-three_prime} and $current->[2] eq q{+} )
    {
        #                       min       max
        # window:                v----+----v                (centered on downstream end)
        # p   |--------------|
        # c                 |---------|
        # n                               |-------------|
        return ( max( $current->[0], $previous->[1], $minimum ),
                 min( $next->[0], $maximum ) );
    }
    else { 
        #             min       max
        # window:      v----+----v                         (centered on upstream end)
        # p                               |-------------|
        # c                 |---------|
        # n   |--------------|
        return ( max( $previous->[1], $minimum ),
                 min( $current->[1], $next->[0], $maximum ) );
    }
}

sub stop_flag_3 {
    croak "Not yet implemented";
}

sub stop_flag_4 {
    croak "Not yet implemented";
}

sub stop_flag_5 {
    croak "Not yet implemented";
}

#==================================================================
# just like stop_flag_2, EXCEPT don't go within -stop_distance of the other end
# (in stop_flag_2, you can go up to the other end).

sub stop_flag_6 {
    my ( $parameters, $previous, $current, $next ) = @_;

    croak "Flag 6 requires the stop distance parameter (-k)"
        unless $parameters->{-stop_distance};

    my ( $minimum, $maximum ) = min_max_distance( $current, $parameters );

    if (   $parameters->{-five_prime} and $current->[2] eq q{-}
        or $parameters->{-three_prime} and $current->[2] eq q{+} )
    {
        # downstream end
        return (
            max($current->[0] + $parameters->{-stop_distance}, $previous->[1],
                $minimum
            ),
            min( $next->[0], $maximum )
        );
    }
    else {
        # upstream end
        return (
            max($previous->[1], $minimum ),
            min($current->[1] - $parameters->{-stop_distance}, $next->[0], $maximum)
        );
    }
}

#==================================================================
# return forward-strand coordinates of the window start/end. don't go into negatives.
# $current = [start,end,strand,attributes]

sub min_max_distance {
    my ( $current, $parameters ) = @_;

    my ( $minimum, $maximum );

    #    start    end
    # 3' <----------- 5'  (5' end on minus strand is 'end'.)
    # 5' -----------> 3'  (3' end on plus strans is also 'end')
    if (   $parameters->{-five_prime} and $current->[2] eq q{-}    # 
        or $parameters->{-three_prime} and $current->[2] eq q{+} ) # 
    {

        $minimum = $current->[1] - $parameters->{-distance};
        $maximum = $current->[1] + $parameters->{-distance};

    }
    #    start    end
    # 5' -----------> 3'  (5' end on plus strand is also 'start')
    # 3' <----------- 5'  (3' end on minus strand is 'start'.)
    else {
        $minimum = $current->[0] - $parameters->{-distance};
        $maximum = $current->[0] + $parameters->{-distance};
    }

    $minimum = ( $minimum < 0 ? 0 : $minimum );

    return ( $minimum, $maximum );
}

__END__


=head1 NAME

 ends_analysis.pl - Produce histogram of GFF data scores, given a GFF annotation

=head1 SYNOPSIS

 ends_analysis.pl -g gene_annotation.gff -b 100 -d 5000 -s 2 -5 -x ID probes.gff

=head1 DESCRIPTION

=head1 OPTIONS

 ends_analysis.pl [OPTION]... [FILE]...

 -g, --gff-annotation GFF 3 annotation file
 -b, --bin-width      histogram bin width                                [100]
 -d, --distance       bp distance from end terminal to search, both ways [5000]
 -s, --stop-flag      when to stop searching (FIXME: explain options)    [2]
 -k, --stop-distance  distance from genes to stop from (flag 6 only)     [1500]
 -3, --three-prime    center analysis on 3' end                          [0]
 -5, --five-prime     center analysis on 5' end                          [1]
 -x, --extract-id     GFF 3 attribute tag pointing to locus ID           [ID]
     --zero           If a bin has no scores mapping to it, but it is 
                      valid for analysis according to stop flag, output 
                      0 instead of 'na'.
     --no-skip        create an entry for gene even if there was nothing
                      found.
 -o, --output         filename to write results to (defaults to STDOUT)
 -v, --verbose        output perl's diagnostic and warning messages
 -q, --quiet          supress perl's diagnostic and warning messages
 -h, --help           print this information
 -m, --manual         print the plain old documentation page

=head2 --stop-flag

=over 

=item 0

For the centered end, go upstream or downstream --distance bases.

=item 2

For the centered end, go upstream or downstream --distance bases but do not
overlap into any other genes.

=item 6

For the centered end, go upstream or downstream --distance bases but do not
overlap OR go within  --stop-distance of into any other genes.

=back

=head1 REVISION

 Version 0.0.1

 $Rev: 396 $:
 $Author: psilva $:
 $Date: 2010-08-13 12:34:16 -0700 (Fri, 13 Aug 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/ends_analysis.pl $:
 $Id: ends_analysis.pl 396 2010-08-13 19:34:16Z psilva $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
