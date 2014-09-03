#!/usr/bin/perl
use warnings;
use strict;
use LWP::UserAgent;
use IPC::Open3;
use File::Copy;
use Getopt::Long;

# perl database_schemas.pl -i /TRIA-NetUtils/rperl change directoryeference_lists/cooke-db-list -o ~/workspace/databases
my ($website_list_file, $apps_dir);
GetOptions(
	'i=s'    => \$website_list_file,
	'a=s'    => \$apps_dir
);
# $website_list_file = "$ENV{\"TRIANetUtils\"}/reference_lists/cooke-lab-site-list" unless defined $website_list_file;
$website_list_file = "/TRIA-NetUtils/reference_lists/cooke-lab-site-list" unless defined $website_list_file;
$apps_dir = '/var/www/apps' unless defined $apps_dir;

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
      }
      else {
	    warn join(": ", $website_name, $response->status_line) . "\n";
	    if($response->status_line =~ m/503 Service Temporarily Unavailable/){
		  my $webapp_dir = join('/', $apps_dir, $website_name);
		  chdir($webapp_dir) or die "Error chdir to $webapp_dir $!";
		  
		  mongrel_stop($website_name);
		  my $restart_status = mongrel_restart($website_name);
		  if($restart_status eq 1){
			my $pid_dir = join('/', $webapp_dir, "tmp", "pids");
			my $pid_files = find_input_files($pid_dir, "pid");

			foreach my $pid_file (sort keys %{$pid_files}){
			      my $file = $pid_files->{$pid_file};
			      warn $file . "\n";
			      copy($file, "$file.old") or die "Copy failed: $!";
			      unlink $file;
			}
			
			chdir($webapp_dir) or die "Error chdir to $webapp_dir $!";;
			mongrel_restart($website_name);
		  }
	    }elsif($response->status_line =~ m/500 Internal Server Error/){
		  my $webapp_dir = join('/', $apps_dir, $website_name);
		  chdir($webapp_dir) or die "Error chdir to $webapp_dir $!";
		  
		  mongrel_stop($website_name);
		  my $restart_status = mongrel_restart($website_name);
		  if($restart_status eq 1){
			my $pid_dir = join('/', $webapp_dir, "tmp", "pids");
			my $pid_files = find_input_files($pid_dir, "pid");

			foreach my $pid_file (sort keys %{$pid_files}){
			      my $file = $pid_files->{$pid_file};
			      warn $file . "\n";
			      copy($file, "$file.old") or die "Copy failed: $!";
			      unlink $file;
			}
			
			chdir($webapp_dir) or die "Error chdir to $webapp_dir $!";;
			mongrel_restart($website_name);
		  }   
	    }
      }
}

sub mongrel_start{	

      my $website_name = shift;
      die "Error lost input file directory" unless defined $website_name;

      warn "Starting server for $website_name....\n";
      my $mongrel_rails_start_cmd  = "mongrel_rails cluster::start";
      warn $mongrel_rails_start_cmd . "\n\n";

      local (*MONGREL_OUT, *MONGREL_IN, *MONGREL_ERR);
      my $pid = open3(\*MONGREL_IN, \*MONGREL_OUT, \*MONGREL_ERR, $mongrel_rails_start_cmd) or die "Error calling open3 for mongrel process start: $!";
      close MONGREL_IN or die "Error closing STDIN to mongrel process start: $!";	

      while (<MONGREL_OUT>){ 
	    chomp $_;
	    warn $_ . "\n";
      }
      close MONGREL_OUT or die "Error closing STDOUT from mongrel process start: $!";

      while (<MONGREL_ERR>){ 
	    chomp $_;
	    warn $_ . "\n";
      }
      close MONGREL_ERR or die "Error closing STERR from mongrel process start: $!";
      wait;
}

sub mongrel_stop{

      my $website_name = shift;
      die "Error lost input file directory" unless defined $website_name;

      warn "Stopping server for $website_name....\n";
      my $mongrel_rails_stop_cmd  = "mongrel_rails cluster::stop";
      warn $mongrel_rails_stop_cmd . "\n\n";

      local (*MONGREL_OUT, *MONGREL_IN, *MONGREL_ERR);
      my $pid = open3(\*MONGREL_IN, \*MONGREL_OUT, \*MONGREL_ERR, $mongrel_rails_stop_cmd) or die "Error calling open3 for mongrel process stop: $!";
      close MONGREL_IN or die "Error closing STDIN to mongrel process stop: $!";	

      while (<MONGREL_OUT>){ 
	    chomp $_;
	    warn $_ . "\n";
      }
      close MONGREL_OUT or die "Error closing STDOUT from mongrel process stop: $!";

      while (<MONGREL_ERR>){ 
	    chomp $_;
	    warn $_ . "\n";
      }
      close MONGREL_ERR or die "Error closing STERR from mongrel process stop: $!";
      wait;
}

sub mongrel_restart{

      my $website_name = shift;
      die "Error lost input file directory" unless defined $website_name;

      warn "Restarting mongrel server for $website_name website app....\n";
      my $mongrel_rails_restart_cmd  = "mongrel_rails cluster::restart";
      warn $mongrel_rails_restart_cmd . "\n\n";

      local (*MONGREL_OUT, *MONGREL_IN, *MONGREL_ERR);
      my $pid = open3(\*MONGREL_IN, \*MONGREL_OUT, \*MONGREL_ERR, $mongrel_rails_restart_cmd) or die "Error calling open3 for mongrel process restart: $!";
      close MONGREL_IN or die "Error closing STDIN to mongrel process restart: $!";	

      while (<MONGREL_OUT>){ 
	    chomp $_;
# 	    warn $_ . "\n";
      }
      close MONGREL_OUT or die "Error closing STDOUT from mongrel process restart: $!";

      my $error_counter = 0;
      while (<MONGREL_ERR>){ 
	    chomp $_;
# 	    warn $_ . "\n";
	    if($_ =~ m/\*\* !!! PID file tmp\/pids\/mongrel.\d+.pid already exists.  Mongrel could be running already.  Check your log\/mongrel.\d+.log for errors./){
		  $error_counter++;
	    }
      }
      close MONGREL_ERR or die "Error closing STERR from mongrel process restart: $!";
      wait;

      my $status = -1;
      if($error_counter > 0){
	    $status = 1;
      }else{
	    $status = 0;
      }

      return $status;
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

