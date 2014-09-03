# serverUtils.pm - 

package serverUtils;
require Exporter;
use strict;
use FileHandle;
use IPC::Open2;
use Term::ReadKey;
use POSIX;
use File::Basename;
use File::Find;
use Cwd qw(chdir abs_path);
use vars qw(@ISA @EXPORT);
@ISA = qw(Exporter);

# Exports all Global Parameters and Commands.
@EXPORT = qw($screen $cpu_cores $VERBOSE startScreenSession checkScreenSession execScreenCmd passwordPrompt getCPUCount

);


# Globals
our ($screen

);

# Edit these paths to ensure that TRIANetUtils can find the executables it needs
$screen			= '/usr/bin/screen';


# Global Parameters
our ($cpu_cores, $VERBOSE

);

# edit these parameters if you want to affect the # cpus TRIANetUtils will use or the assembly parameters
$cpu_cores	= 20;	# -c parameter to newbler
$VERBOSE	= 0;

# These are internal programs within TRIANetUtils and you should not need to edit these unless you are a developer
#$translateStranded		= "$ENV{TRIANetUtils}/bin/translate_stranded.pl";
#$read_ids_file			= 'reads.ids';

# screen related subs
sub startScreenSession {
	my $session_name  = shift or die "lost session name";
	my ($session_id, $session_desc, $isServerOn) = checkScreenSession($session_name);
	if($isServerOn eq "false"){
		warn "Starting new screen session....\n";	

		system($screen,
			'-AmdS',	$session_name
		) == 0 or die "Error: initialization of screen session $session_name failed: $!";
		($session_id, $session_desc, $isServerOn) = checkScreenSession($session_name);
			
		warn $session_desc . "\n";
	}

	return ($session_id, $session_desc, $isServerOn);
}

sub checkScreenSession {
	my $session_name  = shift or die "lost session name";
		
	my ($session_id, $session_desc, $isServerOn);

	$isServerOn = "false";

	my $pid = open2(\*Reader, \*Writer, "$screen -ls");
	close(Writer);
	while(<Reader>){
		chomp $_;
		# 13776.cooke-srv (14-04-03 02:14:49 PM) (Attached|Detached)
		if($_ =~ m/((\d+)\.$session_name\s+\(\d+-\d+-\d+\s+\d+:\d+:\d+\s+(AM|PM)\)\s+\((Attached|Detached)\))/){

			$session_desc = $1;			
			$session_id = $2;
			$session_desc =~ s/\s+/ /g;
			$isServerOn = "true";
		}
	}
	close(Reader);
	
	if($isServerOn eq "false"){
		($session_id, $session_desc, $isServerOn) = ("N/A", "N/A", "false");
	}

	return ($session_id, $session_desc, $isServerOn);
}

sub execScreenCmd{
	my $session_name  = shift or die "lost session name";
	my $session_cmd = shift or die "lost session command";

	my ($session_id, $session_desc, $isServerOn) = checkScreenSession($session_name);
	if($isServerOn eq "true"){
		system($screen,
			'-S', $session_name,
			'-p', 0,
			'-X',
			'stuff', "$session_cmd\n"
		) == 0 or die "Error executing $session_cmd using screen in execScreenCmd: $?";
	}elsif($isServerOn eq "false"){
		warn "Session Name: $session_name does not exist....\n";
		die "Error executing $session_cmd using screen in execScreenCmd: $?";
	}
}

sub passwordPrompt{

	my $password_prompt = shift or die "lost password prompt";
	my $key = 0;
	my $password = "";

	print $password_prompt;

	# Start reading the keys
	ReadMode(3); 

	# Disable the control keys
	while(ord($key = ReadKey(0)) != 10){# This will continue until the Enter key is pressed (decimal value of 10)

		# For all value of ord($key) see http://www.asciitable.com/
		if(ord($key) == 127 || ord($key) == 8){# DEL/Backspace was pressed

			#1. Remove the last char from the password
			chop($password);

			#2 move the cursor back by one, print a blank character, move the cursor back by one
			print "\b \b";

		}elsif(ord($key) < 32){
			# Do nothing with these control characters
		}
		else {
			$password = $password.$key;
 			print "*";
		}
	}

	ReadMode(0); #Reset the terminal once we are done
	#print "\n\nYour super secret password is: $password\n"; 
	return $password;
}


# getCPUCount - Gets the CPU count on the server.

# Output Parameters:                                                        
#                   $CPU_count - Percentage of CPUs on the server.            
sub getCPUCount{
	my ($CPU_count_cmd, $pid, $CPU_count);

	$CPU_count_cmd = "cat /proc/cpuinfo | grep cpu cores";
	
	# getCPUCount I/O command
	local (*NUM_CPU_OUT, *NUM_CPU_IN);
	$pid = open2(\*NUM_CPU_OUT,\*NUM_CPU_IN, $CPU_count_cmd) or die "Error calling open2 in getCPUCount: $!";
	close NUM_CPU_IN or die "Error closing STDIN to getNUM_CPU process: $!";	
	
	print <NUM_CPU_OUT>;
	
	close NUM_CPU_OUT or die "Error closing STDOUT from getNUM_CPU process: $!";
	wait;	
	return $CPU_count;
}

1;
