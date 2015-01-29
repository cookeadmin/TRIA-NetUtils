#!/usr/bin/perl
use warnings;
use strict;
use LWP::UserAgent;
use Getopt::Long;

# perl database_schemas.pl -i /TRIA-NetUtils/rperl change directoryeference_lists/cooke-db-list -o ~/workspace/databases
my ($website_list_file, $apps_dir);
GetOptions(
	'i=s'    => \$website_list_file,
);

# $website_list_file = "$ENV{\"TRIANetUtils\"}/reference_lists/cooke-lab-site-list" unless defined $website_list_file;
$website_list_file = "/TRIA-NetUtils/reference_lists/cooke-lab-site-list" unless defined $website_list_file;

my $website_list = get_website_list($website_list_file);
foreach my $website_name (sort keys %{$website_list}){
#       warn $website_name . "\n";
	my $url = $website_list->{$website_name};

	my $ua = LWP::UserAgent->new;
	$ua->timeout(10);
	$ua->env_proxy;

	my $response = $ua->get($url);

	if ($response->is_success) {
		warn "The $website_name is working so moving on.....\n";
	}else {
		warn join(": ", $website_name, $response->status_line) . "\n";
		if($response->status_line =~ m/500 Can't connect to cookelab:\d+ \(Connection refused\)/){
			my $status = system("service apache2 restart") == 0 or die "Error calling service apache2 restart: $?";
			last;
		}
	}
}

sub get_website_list{
	
    	my $website_list_file = shift or die "lost database list file";
    	
    	my %website_list;
	open(INFILE, "<$website_list_file") or die "Couldn't open file $website_list_file for reading, $!";
	while(<INFILE>){
		chomp $_;
# 		warn "$_\n";
		my ($website_name, $website_link) = split(/\t/, $_);
		$website_list{$website_name} = $website_link;
	}
	close(INFILE) or die "Couldn't close file $website_list_file";
	return \%website_list;
}

sub find_input_files {

	# The directory to parse for input files.
	my $input_dir = shift;
	die "Error lost input file directory" unless defined $input_dir;

	# Input file suffix. Can be anything (.tsv, etc.)
	my $suffix = shift;
	die "Error lost input file suffix" unless defined $suffix;

	if(-d $input_dir){ # Check if $input_dir is a directory.
		opendir(DIR,$input_dir) or die "Error opening $input_dir: $!";
		
		warn "Loading the following files.....\n";

		my %input_files;
		while(my $input_file = readdir(DIR)){
			if ($input_file =~ /\.$suffix$/){
			      warn $input_file . "\n";
			      my $infile = join('/', $input_dir, $input_file);
			      $input_files{$input_file} = $infile;
			}
		}
		closedir(DIR) or die "Error closing $input_dir: $!";

		# Reference to a list of input files.
		return \%input_files;
	}
	else{
		die "Error $input_dir does not exist!\n";
	}
}

