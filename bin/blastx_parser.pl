#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my ($blastx_infile, $blastx_outfile);
GetOptions(
      'i=s'    => \$blastx_infile,
      'o=s'    => \$blastx_outfile,
);

usage() unless (
      defined $blastx_infile
      and defined $blastx_outfile
);

sub usage {

die <<"USAGE";

Usage: $0 -i blastx_infile -o blastx_outfile

Description - 

OPTIONS:

      -i blastx_infile -

      -o blastx_outfile -

USAGE
}

my %blastx_hits = ();
my $i = 0;
open(INFILE,"<$blastx_infile") or die "Error opening $blastx_infile: $!";
while(<INFILE>){
	chomp $_;
	if($i ne 0){
		print $_ . "\n";
		my @split_blastx_entry = split(/\t/, $_);
		my $query_id = $split_blastx_entry[0];
		push(@{$blastx_hits{$query_id}}, $_);
	}
	$i++;
}
close INFILE or die "Error closing $blastx_infile: $!";


warn "Generating blastx tsv file....\n";
open(OUTFILE, ">$blastx_outfile") or die "Couldn't open file $blastx_outfile for writting, $!";
print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch",
      "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
foreach my $query_id (sort keys %blastx_hits){
	for(my $i = 0; $i < 10; $i++){
		print OUTFILE @{$blastx_hits{$query_id}}[$i] . "\n" if(defined(@{$blastx_hits{$query_id}}[$i]));
	}
}
close(OUTFILE) or die "Couldn't close file $blastx_outfile";
