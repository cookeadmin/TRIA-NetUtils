#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Declaration of parameter options to be initialized in @options argument strings.
my $infile;
my @options = (
	'i=s',	\$infile,
);
&GetOptions(@options);

usage() unless (
	defined $infile
);

# Usage details on input parameters for the hierarchial clustering of shared OTUs between multiple libraries.
sub usage {

die <<"USAGE";

Usage: $0 -i infile
OPTIONS:

-i infile

USAGE
}

open INFILE, "<$infile" or die "Error opening $infile for reading: $!";
while(<INFILE>){
      chomp $_;

      warn $_ . "\n";
      system($_) == 0 or die "Error calling $_: $?";
}
close INFILE or die "Error closing $infile: $!";
