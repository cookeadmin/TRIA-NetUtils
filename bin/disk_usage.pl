#!/usr/bin/perl
use warnings;
use strict;

use IPC::Open2;

my $cmd = "df -h";
warn "$cmd\n";

local (*OUTPUT, *INPUT);
my $pid = open2(\*OUTPUT,\*INPUT, $cmd) or die "Error calling open2: $!";
close INPUT or die "Error closing STDIN to df process: $!";	

my $i = 0;
my $header;
my ($size, $used, $free) = 0;
while (<OUTPUT>){
	chomp $_;
	$_ =~ s/\s+/\t/g;
	if($i ne 0){
		if($_ !~ m/mnt|media|^\/\//){
			#print $_ . "\n";
			my @fileSystemInfo = split(/\t/, $_);
			
			$size += $fileSystemInfo[1];
			$used += $fileSystemInfo[2];
			$free += $fileSystemInfo[3];
		
		}
	}
	$i++;
}

print join("\t", "Size", "Used", "Free") . "\n";
print join("\t", $size, $used, $free) . "\n";

close OUTPUT or die "Error closing STDOUT from df process: $!";
wait;
