#!/usr/bin/env perl 

use Modern::Perl '2011';
use autodie;

use Smart::Comments;
use Getopt::Euclid qw( :vars );

use Tie::IxHash;
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
my %best_hit_for;
my %low_hits_for;

for my $file (sort @blast_reports) {
    ### Processing: $file
    
    my $report = Table->new( file => file($file) );
    my $method = 'next_hit';

    my $curr_query;
    my $curr_org;

    HIT:
    while ( my $hit = $report->$method ) {
     
        $curr_query = $hit->query_id;
        $curr_org   = $hit->hit_id;
        
        # skip when hit is myself    
        my ( $query_id, $hit_id ) = map { SeqId->new( full_id => $_ ) } $curr_query, $curr_org;
        next if $query_id->taxon_id eq $hit_id->taxon_id;
     
        my $group = $classifier->classify($curr_org) // 'others';
        
        if ( $best_hit_for{$curr_query}{$curr_org}{$group} ) {
            $low_hits_for{$curr_query}{$curr_org}{$group} = $hit;
            next HIT;
        }
     
        $best_hit_for{$curr_query}{$curr_org}{$group} = $hit;
        next HIT;
    }
    ### Done processing: $file
}
### Done processing all reports

### %best_hit_for
### %low_hits_for

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
