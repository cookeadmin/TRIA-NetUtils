#!/usr/bin/perl
use warnings;
use strict;

use lib "$ENV{'TRIANetUtils'}/lib";	# this allows us to find the module using the lone environment variable we require to be configured by the user
use serverUtils;
# Resource for SSH login without password.
# http://www.thegeekstuff.com/2008/11/3-steps-to-perform-ssh-login-without-password-using-ssh-keygen-ssh-copy-id/

# Resource for Agent admitted failure to sign using the key error message.
# http://www.cyberciti.biz/faq/unix-appleosx-linux-bsd-agent-admitted-failuretosignusingkey/

my $session_name = "cooke-srv";
my $timeout = 1;

my ($session_id, $session_desc, $isServerOn) = serverUtils::checkScreenSession($session_name);

if($isServerOn eq "false"){

	($session_id, $session_desc, $isServerOn) = serverUtils::startScreenSession($session_name);
	
	# screen -S $session_name -p 0 -X stuff 'ssh -X -o TCPKeepAlive=yes -o ServerAliveInterval=99 cookeadmin\@cooke-srv\n'
	serverUtils::execScreenCmd($session_name, "ssh -X -o TCPKeepAlive=yes -o ServerAliveInterval=99 cookeadmin\@cooke-srv");


	# screen -S $session_name -p 0 -X stuff 'serverUtils::passwordPrompt("cookeadmin\@cooke-srv\'s password: ")\n'
	#serverUtils::execScreenCmd($session_name, serverUtils::passwordPrompt("cookeadmin\@cooke-srv\'s password: "));
	sleep($timeout);

	# screen -S $session_name -p 0 -X stuff 'sudo ssh -X -N -o TCPKeepAlive=yes -o ServerAliveInterval=99 -Lcooke-srv.ccis.ualberta.ca:80:cookelab.ccis.ualberta.ca:80 cooke_lab\@cookelab\n'
	serverUtils::execScreenCmd($session_name, "sudo ssh -X -N -o TCPKeepAlive=yes -o ServerAliveInterval=99 -Lcooke-srv.ccis.ualberta.ca:80:cookelab.ccis.ualberta.ca:80 cooke_lab\@cookelab");
	
	# screen -S $session_name -p 0 -X stuff 'serverUtils::passwordPrompt("\n[sudo] password for cookeadmin: ")\n'
	serverUtils::execScreenCmd($session_name, serverUtils::passwordPrompt("\n[sudo] password for cookeadmin: "));
	sleep($timeout);


	# screen -S $session_name -p 0 -X stuff 'serverUtils::passwordPrompt("\ncooke_lab\@cookelab\'s password: ")\n'
	serverUtils::execScreenCmd($session_name, serverUtils::passwordPrompt("\ncooke_lab\@cookelab\'s password: "));
}elsif($isServerOn eq "true"){
	warn "Session Name: $session_name is currently running....\n";
	warn $session_desc . "\n\n";
}

exit 0;
