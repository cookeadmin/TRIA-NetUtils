#!/usr/bin/perl
use warnings;
use strict;
use LWP::UserAgent;

# perl check_CLCServer_status.pl

my $url = "http://snipper.biology.ualberta.ca:7777/CLCServer";

my $ua = LWP::UserAgent->new;
$ua->timeout(10);
$ua->env_proxy;

my $response = $ua->get($url);

if ($response->is_success) {
	warn "The CLCServer is online.....\n";
}else {
	warn join(": ", "CLCServer", $response->status_line) . "\n";
	if($response->status_line =~ m/500 Can't connect to snipper.biology.ualberta.ca:7777 \(Connection refused\)/){
		my $status = system("service CLCGenomicsServer restart") == 0 or die "Error calling service CLCGenomicsServer restart: $?";

	}
}


