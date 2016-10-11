#!/usr/bin/perl -w

use strict;
use lib qw(/share32/perl/);
use Getopt::Long qw (:config gnu_getopt auto_version);
use Bio::SearchIO;

my $P = "/home/db/hapto_blast";

# default values for command-line options
my $g1 = "";
my $g2 = "";
my $delta = 10;
my $minperc = 50;
my $cover = 50;
my $clfile = "";
my $blast = "none";
my $db1 = "$P/blastdb/hyp3v2";
my $db2 = "$P/blastdb/dbali";
my $threshold = 95;

my $parms = GetOptions (
	"g1=s" => \$g1,
	"g2=s" => \$g2,
	"delta=i" => \$delta,
	"minperc=i" => \$minperc,
	"cover=i" => \$cover,
	"clust=s" => \$clfile,
	"blast=s" => \$blast,
	"db1=s" => \$db1,
	"db2=s" => \$db2,
	"threshold=i" => \$threshold
);

# print usage if needed
if (!$parms || @ARGV != 1) {
	print STDERR "Usage: % $0 <file.digest>"
		. " [--g1=<first-group>] [--g2=<second-group>]"
		. " [--delta=<e|%>] [--minperc=<%>] [--cover=<%>]"
		. " [--clust=<file.clust>]"
		. " [--blast=<none|blastx|blastp>"
		. " --db1=<path/db-cds> --db2=<path/db-ali>"
		. " [--threshold=<%>] ]\n";
	exit (1);
}

# check parameter values
if ($blast ne "none" && $blast ne "blastx" && $blast ne "blastp") {
	die ("ERROR! Invalid value for option --blast: $blast!\n");
}
if ($blast ne "none") {
	die ("ERROR! --db1 option mandatory for --blast mode!\n") if (!$db1);
	die ("ERROR! --db2 option mandatory for --blast mode!\n") if (!$db2);
}

# Build exclusion hash
my $exhash = $clfile ? build_exclusion_hash ($clfile) : {};

# get main options from command line
my $infile = shift;

# open infile
open (IN, "<$infile") or
	die ("ERROR! Unable to open digest file '$infile': $!\n");

# build regex from parameters
$g1 = quotemeta ($g1);
$g2 = quotemeta ($g2);
my $size_tok = "\\((\\d+) aa\\)";
my $sort_tok = "\\[([0-9\\.]+)\\]";
my $regex = "^# (\\S+) $size_tok.*?$sort_tok.*?"
	. "with ($g1\\S*) $size_tok $sort_tok.*?"
	. "vs ($g2\\S*) $size_tok $sort_tok";

print STDERR "Parsing digest...\n";
my $count = 0;
my %candidates = ();
my @seq2blast = ();
while (my $line = <IN>) {
	chomp ($line);

	# ensure counting of all not excluded queries
	if ($line =~ /^# --query=(\S+)/) {
		my $query = extract ('\|', $1, 2);
		$count++ if (!defined ($exhash->{$query}));
	}
	
	# get all relevant data from header
	elsif ($line =~ /$regex/) {
		my $query = $1;
		my $length = $2;
		my $delnn = $3;
		my $group1 = $4;
		my $gr1len = $5;
		my $gr1perc = $6;
		my $group2 = $7;
		my $gr2len = $8;
		my $gr2perc = $9;

		# skip duplicate queries (as defined in optional blastclust file)
		next if (defined ($exhash->{$query}));
		
		# compute coverage
		my $len_ratio = 100.0 * $gr1len / $length;
		
		# process only genes passing all three filters
		if ($delnn > $delta && $gr1perc > $minperc && $len_ratio > $cover) {
		
			# extract seq id for first hit (can be from an ignored group)
			<IN>; $line = <IN>;
			chomp ($line);
			my $seqid = extract ('\|', extract ("\t", $line, 6), 2);
			push @seq2blast, $seqid;
			$candidates{$query} = [$length, $delnn,
				$gr1len, $len_ratio, $gr1perc, $group1, $group2, $seqid];
		}
	}
}

my %matches = ();
if ($blast ne "none") {
	
	print STDERR "BLASTing candidate sequences...\n";
	
	# open output file
	my $outfile = "$infile.forblast";
	open (OUT, ">$outfile")
		or die "Unable to write FASTA file '$outfile': $!.\n";

	# fetch seqs from db-cds and write them to FASTA file
	foreach my $id (@seq2blast) {
		my $seq = `fastacmd -d$db1 -s\'$id\'`;
		print OUT $seq;
	}
	close (OUT);
	
	# perform BLASTX vs db-ali
	my $report = "$outfile.$blast";
	`blastall -p$blast -d$db2 -i$outfile -o$report -e1e-10 -FF -b5 -v5`;
	
	# Read BLAST report and loop through queries/results
	my $in = new Bio::SearchIO (-file => $report, -format => 'blast');
	while (my $res = $in->next_result) {
		if (my $hit = $res->next_hit) {
			if (my $hsp = $hit->next_hsp) {
				
				# store hit if identity > threshold
				if ($hsp->frac_identical * 100.0 > $threshold) {
					my $seqid = extract ('\|', $res->query_name, 2);
					$matches{$seqid} = $hit->name;
				}
			}
		}
	}
}

print STDERR "Printing results...\n";
printf "# $infile [%d proteins retained out of %d] --delta=%.2f --minperc=%d"
	. " --cover=%d", scalar keys %candidates, $count, $delta, $minperc, $cover;
printf " --blast=$blast --threshold=%2.f", $threshold if ($blast ne "none");
printf " --clust=%s", $clfile if ($clfile);
print "\n";
print "# query_id\tlength\tdelta\tlen_1\tcover_1\tperc_1\tgroup_1         "
	. "\tgroup_2         \t1st_hit   ";
print "\tali_file" if ($blast ne "none");
print "\n";
foreach my $query (sort keys %candidates) {
	my ($length, $delnn, $gr1len, $len_ratio,
		$gr1perc, $group1, $group2, $seqid) = @{$candidates{$query}};
	printf "%-10s\t%d\t%.2f\t%d\t%.2f\t%.2f\t%-16s\t%-16s\t%-10s", $query, $length,
		$delnn, $gr1len, $len_ratio, $gr1perc, $group1, $group2, $seqid;
	if (defined ($matches{$seqid})) {
		print "\t" . extract ('\|', $matches{$seqid}, 1);
	}
	print "\n";
}

close (IN);



# parse blastclust output and build exclusion hash
sub build_exclusion_hash {
	my ($list) = @_;
	my %hash = ();
	open (IN, "<$list") or
		die ("ERROR! Unable to open blastclust file '$list': $!\n");
		
	# process each cluster (one per line)
	while (my $line = <IN>) {

		# extract ids from current cluster
		chomp ($line);
		my @ids = split (" ", $line);
		
		# discard first id (will be used as unique copy)
		shift @ids;
		
		# store remaining ids into exclusion hash
		while (my $id = shift @ids) {
			$hash{$id} = 1;
		}
	}
	close (IN);

	return \%hash;
}



sub extract {
	my ($delim, $str, $field) = @_;

	my @subs = split ($delim, $str);
	return $subs[$field];
}
