#!/usr/bin/perl -w

use strict;
use lib qw(/share32/perl/);
use Getopt::Long qw (:config gnu_getopt auto_version);

my $P = "/home/db/hapto_blast";

# default values for command-line options
my $idfile = "$P/data/codesv3.index";
my $grfile = "";
my $minlen = 100;
my $exclude = "self";
my $sort = "perc";
my $style = "short";
my $desc = 0;
my $db = "$P/blastdb/hyp3v2";

my $parms = GetOptions (
	"ids=s" => \$idfile,
	"groups=s" => \$grfile,
	"minlen=i" => \$minlen, 
	"exclude=s" => \$exclude,
	"sort=s" => \$sort,
	"style=s" => \$style,
	"desc" => \$desc,
	"db=s" => \$db
);

# Print usage if needed
if (!$parms || @ARGV != 1) {
	print STDERR "Usage: % $0 <report.m9>"
		. " --groups=<groups.otu> [--ids=<orgs.index>]"
		. " [--minlen=<aa>] [--exclude=<none|self|group>]"
		. " [--sort=<perc|eval>] [--style=<short|long>]"
		. " [--desc --db=<path/db-ali>]\n";
	exit (1);
}

# Check parameter values
if (!$grfile) {
	die ("ERROR! --groups option mandatory!\n");
}
if ($exclude ne "none" && $exclude ne "self" && $exclude ne "group") {
	die ("ERROR! Invalid value for option --exclude: $exclude!\n");
}
my $sortf;
if ($sort eq "eval") {
	$sortf = 0;
} elsif ($sort eq "perc") {
	$sortf = 1;
} else {
	die ("ERROR! Invalid value for option --sort: $sort!\n");
}
if ($style ne "short" && $style ne "long") {
	die ("ERROR! Invalid value for option --style: $style!\n");
}
if ($desc && !$db) {
	die ("ERROR! --db option mandatory for --desc mode!\n");
}

# Build organism and group lookups
my $org_ids = build_lookup ($idfile, "Unkn", "Unknown organism");
my $grp_ids = build_otu_lookup ($grfile, "Unkn", "Unknown group");

# Get filenames from command line
my $infile = shift;

open (IN, "<$infile") or
	die ("ERROR! Unable to open BLAST table file '$infile': $!\n");

# Parse BLAST table
my %top_res = ();
my %low_res = ();
my $query_name = "";
my $query_length = "";
my $query_org = "";
my $query_grp = "";
while (my $line = <IN>) {

	# skip comments and empty lines
	next if ($line =~ /^#/ || $line =~ /^\s*$/);

	# get relevant fields for HSP
	chomp ($line);
	my @fields = split ('\t', $line);
	my $query_id = $fields[0];
	my $hit_id = $fields[1];
	my $percent = $fields[2];
	my $length = $fields[3];
	my $evalue = $fields[10];

	# get query length and id from first hit
	$query_length = $length if (!$query_length);
	$query_name = $query_id if (!$query_name);
	
	# get query name and group using indexes
	if (!$query_org) {
		my $query_code = substr (extract ('\|', $query_id, 2), 0, 4);
		$query_org = $org_ids->{$query_code};
		die ("ERROR! $query_code not defined in '$idfile'!\n")
			if (!defined ($query_org));
		$query_grp = $grp_ids->{$query_org};
		die ("ERROR! $query_org not defined in '$grfile'!\n")
			if (!defined ($query_grp));		
	}

	# skip short queries and short HSPs (absolute min-len)
	last if ($query_length < $minlen);
	next if ($length < $minlen);
	
	# get organism id from hit id
	$hit_id =~ m/\S+\|\S+\|([A-Za-z]+)\d+/;
	my $org_id = $1;
	
	# get organism name and group using indexes
	my $org = $org_ids->{$org_id};
	die ("ERROR! $org_id not defined in '$idfile'!\n") if (!defined ($org));
	my $group = $grp_ids->{$org};
	die ("ERROR! $org not defined in '$grfile'!\n") if (!defined ($group));
	
	# optionally exclude query org and/or group from ranking
	my $exclf = 0;
	if ($exclude eq "self") {
		$exclf = 1 if ($org eq $query_org);
	} elsif ($exclude eq "group") {
		$exclf = 1 if ($group eq $query_grp);
	}
	$group = "!" . $group if ($exclf && substr ($group, 0, 1) ne "!");
		
	# fetch sequence description if needed
	if ($desc) {
		`fastacmd -d $db -s "$hit_id"` =~ /^>\S+\s+\[.*?\]\s+(.*)\n/;
		$hit_id = sprintf ("%-20s\t%s", $hit_id, $1);
	}
	
	# store HSP details (extract first organism hits from remaining hits)
	if (!defined ($top_res{$org})) {
		$top_res{$org} = [$evalue, $percent, $length, $hit_id, $group];
	} elsif (!defined ($low_res{$hit_id})) {		# store only first HSP
		$low_res{$hit_id} = [$evalue, $percent, $length, $org, $group];
	}
}
close (IN);

# Exit if no hit (this will produce no output)
exit if (!$query_name);

# Build header
my $head = "# --query=$query_name --groups=$grfile --minlen=$minlen"
	. " --exclude=$exclude --sort=$sort\n";

# Sort top results according to evalue or percent identity
my @tops = sort { ${$top_res{$a}}[$sortf] <=> ${$top_res{$b}}[$sortf] }
	keys %top_res; @tops = reverse @tops if ($sortf);

# Build top part of result table...
my $format = "%4s\t%-16s\t%-8s\t%-8s\t%-8s\t%-32s\t%s\n";
my %grp_hash = ();
my @grp_array = ();
my $body = "";
foreach my $org (@tops) {
	my ($evalue, $percent, $length, $hit_id, $group) = @{$top_res{$org}};

	# ... and store best evalue/percent for each group
	my $rank = "";		# discard groups starting with "!" char
	if (substr ($group, 0, 1) ne "!" and !defined ($grp_hash{$group})) {
		$grp_hash{$group} = [$evalue, $percent, $length];
		push @grp_array, $group;
		$rank = @grp_array;
	}
	
	$body .= sprintf $format,
		$rank, $group, $evalue, $percent, $length, $org, $hit_id;
	
}

# Compute delta between 1st and 2nd group
if (@grp_array > 1) {
	my $a = ${$grp_hash{$grp_array[0]}}[$sortf];
	my $b = ${$grp_hash{$grp_array[1]}}[$sortf];
	if (!$sortf) {	# evalue mode
		$a = $a > 0 ? (- log ($a) / log (10)) : 300;
		$b = $b > 0 ? (- log ($b) / log (10)) : 300;
	}
	my $delnn = $a - $b;

	my @id_fields = split ('\|', $query_name);	# get basename
	$head .= sprintf (
		"# %s (%d aa) 1st-to-2nd-delta [%.2f]"
		. " with %s (%d aa) [%.2f] vs %s (%d aa) [%.2f]\n",
		$id_fields[2], $query_length, $delnn,
		$grp_array[0], ${$grp_hash{$grp_array[0]}}[2], $a,
		$grp_array[1], ${$grp_hash{$grp_array[1]}}[2], $b
	);
}

# Build remaining of header
$head .= sprintf "# %-2s\t%-16s\t%-8s\t%-8s\t%-8s\t%-32s\t%-20s%s\n",
	"rank", "group", "evalue", "identity", "length", "organism", "accession",
	($desc ? "\tdescription" : "");

# Output head and top part of results
print $head;
print $body;
printf "# remaining hits (%d)\n", scalar keys %low_res;

# Append optional lower part of results
exit if ($style eq "short");

# Sort lower results according to evalue or percent identity
my @lows = sort { ${$low_res{$a}}[$sortf] <=> ${$low_res{$b}}[$sortf] }
	keys %low_res; @lows = reverse @lows if ($sortf);
	# Note: potential gotcha when several groups have an e-value of 0
	# in such cases, the querying group may not get the first rank
	# which can disrupt the downstream filtering process

# Build lower part of result table
foreach my $hit_id (@lows) {
	my ($evalue, $percent, $length, $org, $group) =
		@{$low_res{$hit_id}};
	printf $format,
		"", $group, $evalue, $percent, $length, $org, $hit_id;
}



# parse indexes and build look-up hash
# indexes can be either id->id or id->field or id->field1/field2
# unknown parameters specify the default key->value pair
sub build_lookup {
	my ($index, $key, $value) = @_;
	my %lookup = ($key => $value);
	open (IN, "<$index") or
		die ("ERROR! Unable to open index file '$index': $!\n");
	while (my $line = <IN>) {
		chomp ($line);
		my ($id, $field1, $field2) = split ('\t', $line);
		if (defined ($field2)) {
			$lookup{$id} = [$field1, $field2];
		} else {
			$lookup{$id} = $field1;
		}
	}
	close (IN);

	return \%lookup;
}


sub build_otu_lookup {
	my ($index, $key, $value) = @_;
	my %lookup = ();
	open (IN, "<$grfile") or
		die ("ERROR! Unable to open otu file '$grfile': $!\n");
	while (my $line = <IN>) {
		# skip comments and empty lines
		next if ($line =~ /^#/ || $line =~ /^\s*$/);
		
		# get group name and corresponding organisms
		chomp ($line);
		my ($grp, $org_str) = split (':', $line);
		$grp =~ s/^\s+|\s+$//g;			# trim spaces
		my @tops = split (',', $org_str);
		
		# store organisms
		foreach my $org (@tops) {
			$org =~ s/^\s+|\s+$//g;		# trim spaces
			$lookup{$org} = $grp;
		}
	}
	close (IN);
	
	return \%lookup;
}


sub extract {
	my ($delim, $str, $field) = @_;

	my @subs = split ($delim, $str);
	return $subs[$field];
}
