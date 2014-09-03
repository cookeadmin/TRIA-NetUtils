#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use POSIX;
# sudo perl pg_dump_databases.pl -i /TRIA-NetUtils/reference_lists/cooke-db-list -o /home/cookeadmin/workspace

my ($db_file, $output_dir);
GetOptions(
	'i=s'    => \$db_file,
	'o=s'    => \$output_dir
);

usage() unless (
      defined $output_dir
);

$db_file = "$ENV{\"TRIANetUtils\"}/reference_lists/cooke-db-list" unless defined $db_file;

my ($pg_dump, $tar);
$pg_dump	= '/usr/bin/pg_dump';
$tar		= '/bin/tar';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i db_file -o output_dir
    
Description - 
    
OPTIONS:
     
      -i db_file - 
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my ($db_super_user, $db_password, $date_stamp);

$db_super_user = "postgres";
$db_password = "postgresql*ygg";

$date_stamp = strftime("%Y-%m-%d", localtime);

my $db_dirname = join("-", "postgresdb", $date_stamp);

my $db_output_dir = join('/', $output_dir, $db_dirname);

# Create output directory if it doesn't already exist.
unless(-d $db_output_dir){
	mkdir($db_output_dir, 0777) or die "Can't make directory: $!";
}

$ENV{'PGPASSWORD'} = $db_password;

warn $db_file . "\n";
my $db_name_list = get_db_names($db_file);
my %sqlfiles;
foreach my $db_name (@{$db_name_list}){
	warn "Dumping database $db_name to $db_output_dir\n";

	my $sqlfilename = join("_", $db_name, $date_stamp) . ".psql";
	my $sqlfile = join('/', $db_output_dir, $sqlfilename);
	warn $sqlfile . "\n\n";

	system($pg_dump,
		'--clean',
		'-U', $db_super_user,
		$db_name,
		'-w', 
		'-f', $sqlfile
	) == 0 or die "Error executing pg_dump: $?";

	$sqlfiles{$sqlfilename} = $sqlfile;
}



warn "Compressing $db_dirname using tar\n\n";
chdir($output_dir) or die "Cant chdir to $output_dir $!";
my $tarball_outfile = join('/', $output_dir, "$db_dirname.tar.gz");
system($tar,
	'-zcvf', 
	$tarball_outfile,
	$db_dirname
) == 0 or die "Error executing tar -zcvf $tarball_outfile $db_dirname: $?";


sub get_db_names{
	
    	my $db_list = shift or die "lost database list file";
    	
    	my @db_names;
	open(INFILE, "<$db_list") or die "Couldn't open file $db_list for reading, $!";
	while(<INFILE>){
		chomp $_;
		#warn "$_\n";
		push(@db_names, $_);
	}
	close(INFILE) or die "Couldn't close file $db_list";
	return \@db_names;
}

exit 0;
