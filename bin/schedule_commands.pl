#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use lib "$ENV{'TRIANetUtils'}/lib";	# this allows us to find the module using the lone environment variable we require to be configured by the user
use serverUtils;

my $cmd_list_file;
GetOptions(
	'i=s'    => \$cmd_list_file
);
$cmd_list_file = "$ENV{\"TRIANetUtils\"}/reference_lists/screen-commands-list" unless defined $cmd_list_file;

my $command_list = get_command_list($cmd_list_file);

foreach my $session_name (keys %{$command_list}){
	warn $command_list->{$session_name} . "\n";
	my ($session_id, $session_desc, $isServerOn) = serverUtils::checkScreenSession($session_name);

	if($isServerOn eq "false"){

		($session_id, $session_desc, $isServerOn) = serverUtils::startScreenSession($session_name);
	
		serverUtils::execScreenCmd($session_name, $command_list->{$session_name});

	}elsif($isServerOn eq "true"){
		warn "Session Name: $session_name is currently running....\n";
		warn $session_desc . "\n\n";
	}
}



sub get_command_list{
	
    	my $cmd_list_file = shift or die "lost command list file";
    	
    	my %command_list;
	open(INFILE, "<$cmd_list_file") or die "Couldn't open file $cmd_list_file for reading, $!";
	while(<INFILE>){
		chomp $_;
		#warn "$_\n";
		my @split_entry = split(/\t/, $_);
		my ($session_name, $cmd) = ($split_entry[0], $split_entry[1]);
		#warn join("\n", $session_name, $cmd);
		$command_list{$session_name} = $cmd;

	}
	close(INFILE) or die "Couldn't close file $cmd_list_file";
	return \%command_list;
}

exit 0;
