#!/usr/bin/env perl 

use Modern::Perl '2011';
use autodie;

use Smart::Comments;
use Getopt::Euclid qw( :vars );

use Tie::IxHash::Easy;
use File::Basename;
use File::Find::Rule;
use Path::Class 'file', 'dir';
use List::AllUtils 'max', 'firstidx';

use Bio::MUST::Core;
use aliased 'Bio::MUST::Core::Ali';
use aliased 'Bio::MUST::Core::SeqId';
use aliased 'Bio::MUST::Core::IdList';
use aliased 'Bio::MUST::Core::IdMapper';
use Bio::MUST::Core::Utils 'change_suffix';
use aliased 'Bio::FastParsers::Blast::Table';
use aliased 'Bio::MUST::Core::Taxonomy';
use aliased 'Bio::MUST::Core::Taxonomy::Filter';
use aliased 'Bio::MUST::Core::Taxonomy::Criterion';
use aliased 'Bio::MUST::Core::Taxonomy::Category';
use aliased 'Bio::MUST::Core::Taxonomy::Classifier';



# build taxonomy objects
my $tax = Taxonomy->new( tax_dir => $ARGV_taxdir );

### Processing OTUs...
open my $in, '<', file($ARGV_otu_file);

# build classifier from labels
my @categories;

while ( my $line = <$in> ) {    
    chomp $line;
    ### $line
    my ($label, $otu) = split ':', $line;
    my $list = IdList->new( ids => [ split ',', $otu ] );
    ### $label
    ### $otu
    ### $list
    my $filter    = $tax->tax_filter( $list );
    my $criterion = Criterion->new( tax_filter => $filter );
    my $category  = Category->new(
        label    => $label,
        criteria => [ $criterion ],
    );
    push @categories, $category;
}

my $classifier = Classifier->new( categories => \@categories );
### $classifier



### Processing: $ARGV_reports

my $blast_suffix = qr{\. $ARGV_suffix \z}xmsi;

my @blast_reports = File::Find::Rule
    ->file()
    ->name($blast_suffix)
    ->maxdepth(1) 
    ->in($ARGV_reports)
;
### @blast_reports

### Processing blast reports...
for my $file (sort @blast_reports) {
    ### Processing: $file

    my $report = Table->new( file => file($file) );
    my $method = 'next_hit';
    
    my %best_for;
    tie my %best_hit_for, 'Tie::IxHash::Easy';
    tie my %low_hits_for, 'Tie::IxHash::Easy';
#    my $curr_query;
#    my $curr_org;

    # Parse table and store best hit for a query_org/hit_org pair in a hash
    # and the remaining low scoring hits in another hash
    HIT:
    while ( my $hit = $report->$method ) {

        my $query_id = $hit->query_id;
        my $hit_id   = $hit->hit_id;
        ### $query_id
        ### $hit_id
        ### eval: $hit->evalue
        
        # skip when hit is myself    
        my ( $query_taxid, $hit_taxid ) = map { $_->taxon_id } 
                                          map { SeqId->new( full_id => $_ ) } $query_id, $hit_id
                                          ;
        next if $query_taxid eq $hit_taxid;
     
        my $group = $classifier->classify($hit_id) // 'others';
        
        unless ( $best_for{$query_id}{$hit_taxid} ) {

            $best_hit_for{$query_id}{$hit_taxid}{hit_id} = $hit_id;
            $best_hit_for{$query_id}{$hit_taxid}{group}  = $group;
            $best_hit_for{$query_id}{$hit_taxid}{hit}    = $hit;

            $best_for{$query_id}{$hit_taxid} = 1;
            next HIT;
        }

        $low_hits_for{$query_id}{$hit_taxid}{$hit_id}{group} = $group;
        $low_hits_for{$query_id}{$hit_taxid}{$hit_id}{hit}   = $hit;
        next HIT;
    }

    my $basefile      = basename($file);
    my $outfile_best  = change_suffix($basefile, '.best');
    my $outfile_lows  = change_suffix($basefile, '.lows');
    my $outfile_delta = change_suffix($basefile, '.delta');
    ### $outfile_best 
    
    # Write table for the best scoring hits
    open my $out_best, '>', $outfile_best;
    for my $query_id (keys %best_hit_for) {
        for my $hit_taxid (keys %{ $best_hit_for{$query_id} }) {

            my $hit_id = $best_hit_for{$query_id}{$hit_taxid}{hit_id};
            my $group  = $best_hit_for{$query_id}{$hit_taxid}{group};
            my $hit    = $best_hit_for{$query_id}{$hit_taxid}{hit};

	        my @hit_values = map { $hit-> $_ } qw(evalue percent_identity hsp_length);

            say {$out_best} join "\t", $query_id, @hit_values, $hit_id, $group;
        }
    }
    close $out_best;

    # Write table for low scoring hits
    open my $out_lows, '>', $outfile_lows;
    for my $query_id (keys %low_hits_for) {
        for my $hit_taxid (keys %{ $low_hits_for{$query_id} }) {
            for my $hit_id (keys %{ $low_hits_for{$query_id}{$hit_taxid} }) {
             
                my $group  = $low_hits_for{$query_id}{$hit_taxid}{$hit_id}{group};
                my $hit    = $low_hits_for{$query_id}{$hit_taxid}{$hit_id}{hit}; 

	            my @hit_values = map { $hit-> $_ } qw(evalue percent_identity hsp_length);

                say {$out_lows} join "\t", $query_id, @hit_values, $hit_id, $group;
            }
        }
    }
    close $out_lows;

    # Write table for the best scoring hits
    open my $out_delta, '>', $outfile_delta;
    for my $query_id (keys %best_hit_for) {
        HIT_ID:
        for my $hit_taxid (keys %{ $best_hit_for{$query_id} }) {

            my $group  = $best_hit_for{$query_id}{$hit_taxid}{group};
            next HIT_ID if $group =~ m/^!/xmsg;

            my $hit_id = $best_hit_for{$query_id}{$hit_taxid}{hit_id};
            my $hit    = $best_hit_for{$query_id}{$hit_taxid}{hit};

	        my @hit_values = map { $hit-> $_ } qw(evalue percent_identity hsp_length);

            say {$out_delta} join "\t", $query_id, @hit_values, $hit_id, $group;
        }
    }
    close $out_delta;

#    ### %best_hit_for
#    ### %low_hits_for
    ### Done processing: $file
}
### Done processing all reports

sub comp_delta { 
   my ($a, $b) = @_;
   return $a-$b; 
}

=head1 NAME

new_digest.pl

=head1 VERSION

This documentation refers to new_digest.pl version 0.0.1

=head1 USAGE

new_digest.pl --report-dir=<dir> --otu[-file]=<file> --taxdir=<dir>

=head1 REQUIRED ARGUMENTS

=over

=item --reports=<dir>

Path to blast reports.

=for Euclid:
    dir.type: string

=item --suffix=<str>

BLAST reports file extension.

=for Euclid:
    str.type: string

=item --otu[-file]=<file>

Path to artificial groups' file (user defined).

=for Euclid:
    file.type: string

=back

=head1 OPTIONAL ARGUMENTS

=over

=item --taxdir=<dir>

Path to local NCBI taxonomy DB.

=for Euclid:
    dir.type: string

=item --filesdir=<dir>

Path to ali files directory.

=for Euclid:
    dir.type: string

=item --namesfile=<file>

Path to ali files directory.

=for Euclid:
 
=item --gca-mapper=<file>

Path to ali files directory.

=for Euclid:
    file.type: string

=item --user=<num>

User identifier. 

    00, Mick VV
    01, Denis B
    02, Richard G
    03, Léonard R
    04, Di Franco A
    ...
    
=for Euclid:
    num.type: int

=item --count=<num>

Last defined user taxid. Number from which to start count for new user taxids.

=for Euclid:
    num.type: int

=back
