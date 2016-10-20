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
use Statistics::Descriptive;

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

    # Parse table and store best hit for a query_org/hit_org pair in a hash
    # and the remaining low scoring hits in another hash
    HIT:
    while ( my $hit = $report->$method ) {

        next HIT if $hit->hsp_length < $ARGV_hsp_length;
        next HIT if $hit->evalue     > $ARGV_evalue;

        my $query_id = $hit->query_id;
        my $hit_id   = $hit->hit_id;
        
        # skip when hit is myself    
        my ( $query_taxid, $hit_taxid ) = map { $_->taxon_id } 
                                          map { SeqId->new( full_id => $_ ) } $query_id, $hit_id
                                          ;
        next if $query_taxid eq $hit_taxid;
     
        # classify org accoording to config file
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
    say {$out_best} join "\t", '#query_id', 'evalue', 'percent_identity', 'bit_score', 'hsp_length', 'hit_id', 'group';
    for my $query_id (keys %best_hit_for) {
        for my $hit_taxid (keys %{ $best_hit_for{$query_id} }) {

            my $hit_id = $best_hit_for{$query_id}{$hit_taxid}{hit_id};
            my $group  = $best_hit_for{$query_id}{$hit_taxid}{group};
            my $hit    = $best_hit_for{$query_id}{$hit_taxid}{hit};
        
	        my @hit_values = map { $hit-> $_ } qw(evalue percent_identity bit_score hsp_length);
            ### @hit_values

            say {$out_best} join "\t", $query_id, @hit_values, $hit_id, $group;
        }
    }
    close $out_best;

    # Write table for low scoring hits
    open my $out_lows, '>', $outfile_lows;
    say {$out_lows} join "\t", '#query_id', 'evalue', 'percent_identity', 'bit_score', 'hsp_length', 'hit_id', 'group';
    for my $query_id (keys %low_hits_for) {
        for my $hit_taxid (keys %{ $low_hits_for{$query_id} }) {
            for my $hit_id (keys %{ $low_hits_for{$query_id}{$hit_taxid} }) {
             
                my $group  = $low_hits_for{$query_id}{$hit_taxid}{$hit_id}{group};
                my $hit    = $low_hits_for{$query_id}{$hit_taxid}{$hit_id}{hit}; 

	            my @hit_values = map { $hit-> $_ } qw(evalue percent_identity bit_score hsp_length);

                say {$out_lows} join "\t", $query_id, @hit_values, $hit_id, $group;
            }
        }
    }
    close $out_lows;

    # Write table for the computed deltas
    open my $out_delta, '>', $outfile_delta;
    say {$out_delta} join "\t", '#query_id', 'delta-' . $ARGV_delta_mode . "-$ARGV_score", 'group_1', 'group_2';

    QUERY_ID:
    for my $query_id (keys %best_hit_for) {

        my @groups;
        my @values_a;
        my @values_b;

        HIT_ID:
        for my $hit_taxid (keys %{ $best_hit_for{$query_id} }) {

            my $group  = $best_hit_for{$query_id}{$hit_taxid}{group};
            
            # skip unconsidered groups
            next HIT_ID if $group =~ m/^!/xmsg;
            # consider only 2 groups to compute delta
            last HIT_ID if scalar @groups > 2;

            push @groups, $group unless grep { $_ eq $group } @groups;    
#            ### @groups
            
            my $hit_id = $best_hit_for{$query_id}{$hit_taxid}{hit_id};
            my $hit    = $best_hit_for{$query_id}{$hit_taxid}{hit};
#            ### $hit

	        my ($value, $hsp_len) = map { $hit-> $_ } $ARGV_score, qw(hsp_length);
            push @values_a, $value if scalar @groups == 1;
            push @values_b, $value if scalar @groups == 2;
        }

        my $delta = _comp_delta(\@values_a, \@values_b, $ARGV_delta_mode);
        ### $delta

        say {$out_delta} join "\t", $query_id, $delta, $groups[0], $groups[1] // 'NA'; 
    }
    close $out_delta;

#    ### %best_hit_for
#    ### %low_hits_for
    ### Done processing: $file
}
### Done processing all reports

sub _comp_delta { 
    my ($a, $b, $mode) = @_;

    if ($mode eq 'med') {
        my $stat_a = Statistics::Descriptive::Full->new();
        $stat_a->add_data(@$a); 
        my $stat_b = Statistics::Descriptive::Full->new();
        $stat_b->add_data(@$b); 

        my $med_a = $stat_a->median();
        my $med_b = $stat_b->median() // 0;
        ### $med_a
        ### $med_b

        return $med_a - $med_b
    }
    else {
#        ### a: @$a[0]
#        ### b: @$b[0]
        return @$a[0]-@$b[0]; 
    }
}

=head1 NAME

new_digest.pl

=head1 VERSION

This documentation refers to new_digest.pl version 0.0.1

=head1 USAGE

new_digest.pl --report-dir=<dir> --otu[-file]=<file> --taxdir=<dir>

=head1 REQUIRED ARGUMENTS

=over

=item --taxdir=<dir>

Path to local NCBI taxonomy DB.

=for Euclid:
    dir.type: string

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

=item --score=<str>

Choose either percent_identity or bit_score to compute delta

=for Euclid:
    str.type: str

=back

=head1 OPTIONAL ARGUMENTS

=over

=item --delta[-mode]=<str>

Use top hit or median to compute delta.

=for Euclid:
    str.type: str
    str.default: qw(top)

=item --hsp[-length]=<int>

Filter hits according to a min hsp length.

=for Euclid:
    int.type: int
    int.default: 0

=item --evalue=<num>

Filter hits according to a threshold e-value.

=for Euclid:
    num.type: number
    num.default: 1e-10

    
=back
