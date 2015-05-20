#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;
use DBD::Pg;
use File::Basename;
use POSIX;
use open qw/:std :utf8/;

#  sudo perl lims_database_bulk_upload.pl -i ~/workspace/bulk_upload_files/reformated_csv_files/M004_fungal_associates.csv -n triad_demo -o ~/workspace/triad_lims_demo_backups
my ($input_dir, $webapp_name, $backup_dir);
GetOptions(
	'i=s'    => \$input_dir,
	'n=s'    => \$webapp_name,
	'o=s'    => \$backup_dir
);

usage() unless (
      defined $input_dir
#       and defined $webapp_name
#       and defined $backup_dir
);

# my ($db_name, $webapp_dir);
# if($webapp_name eq 'TRIA-LIMS'){
# 
#       $db_name = 'triad_production';
#       $webapp_dir = '/var/www/apps/triad';
# 
# }elsif($webapp_name eq 'triad-demo'){
# 
#       $db_name = 'tria_lims_production';
#       $webapp_dir = '/var/www/apps/triad_demo';
# 
# }else{
#       warn "\nYou entered $webapp_name. Please re-enter the -n webapp_name parameter (must be triad-demo or TRIA-LIMS).\n";
#       usage();
# }
# 
# my $pg_dump	= '/usr/bin/pg_dump';
# 
# sub usage {
#     
# die <<"USAGE";
#     
# Usage: $0 -i infile -n webapp_name -o backup_dir
#     
# Description - 
#     
# OPTIONS:
#       -i infile - 
#       -n webapp_name - 
#       -o backup_dir -
#     
# USAGE
# }
# 
# # Create output directory if it doesn't already exist.
# unless(-d $backup_dir){
# 	mkdir($backup_dir, 0777) or die "Can't make directory: $!";
# }

# my ($db_host, $db_super_user, $db_password);
# $db_host = "localhost";
# $db_super_user = "postgres";
# $db_password = "postgresql*ygg";

# warn "Preparing some_file for bulk upload into $db_name......" . "\n\n";

# warn "Dumping database $db_name to $backup_dir......\n";

# my $date_stamp = strftime("%Y-%m-%d_%H-%M-%S", localtime);
# my $date_stamp = strftime("%Y-%m-%d", localtime);

# $ENV{'PGPASSWORD'} = $db_password;
# 
# my $sqlfilename = join("_", $db_name, $date_stamp) . ".psql";
# my $sqlfile = join('/', $backup_dir, $sqlfilename);
# warn $sqlfile . "\n\n";

# system($pg_dump,
#       '--clean',
#       '-U', $db_super_user,
#       $db_name,
#       '-w', 
#       '-f', $sqlfile
# ) == 0 or die "Error executing pg_dump: $?";

# maintenance_start($webapp_name, $webapp_dir);
# 
# warn "Connecting to database $db_name\n";

# initialize database handlers.
# my($dbh, $sth);
# 
# $dbh = DBI->connect(
#       "dbi:Pg:dbname=$db_name;host=$db_host", 
#       $db_super_user, 
#       $db_password
# ) or die "Unable to connect: $DBI::errstr\n";

my ($database_files, $database_file_count) = find_files($input_dir, "txt");

foreach my $file_name (sort keys %{$database_files}){
	warn "Processing " . $file_name . ".....\n";
	my $database_files_infile = $database_files->{$file_name};
	open(INFILE, "<$database_files_infile") or die "Couldn't open file $database_files_infile for reading, $!";
	my $i = 0;
	
	
	my $database_tablename = fileparse($file_name, qr/\.txt/);
	my @split_header = ();
	while(<INFILE>){
		chomp $_;
		if($i eq 0){
			@split_header = split(/\t/, $_);
		}
		else{
			my @split_entry = split(/\t/, $_);
			for(my $j = 1; $j < scalar(@split_header); $j++){
				warn join("=", $split_header[$j], "'$split_entry[$j]'") . "\n";
			}
			die;
			warn "UPDATE $database_tablename SET WHERE\n";
		}
		$i++;
	}
	close(INFILE) or die "Couldn't close file $database_files_infile";
# 	
# 	my @table_attributes = ();
# 	foreach my $column_name (keys %column_values){
# 		if($column_name !~ m/^id$/){
#  			push(@table_attributes, $column_name);
# 		}
# 	}
# 	for(my $i = 0; $i < $row_counter; $i++){
# 	
# 	}
# 	warn "UPDATE $database_tablename SET WHERE\n";
# 	die;
# 	UPDATE Customers
# SET ContactName='Alfred Schmidt', City='Hamburg'
# WHERE CustomerName='Alfreds Futterkiste'; 
}

# for(my $i = 0; $i < $row_counter; $i++){
# 
# }


# warn "Disconnecting from database $db_name\n";
# my $rc = $dbh->disconnect() or die "Unable to disconnect: $DBI::errstr\n";
# maintenance_stop($webapp_name, $webapp_dir);

# (\%files, $file_counter) = find_files($infile_dir) - Find all files in the specified input file directory with the file extension *.suffix.
#
# Input paramater(s):
#
# $infile_dir - The input file directory.
#
# $suffix - The file extension suffix.
#
# Output paramater(s):
#
# \%files - A hash reference containing all the files with file extension *.suffix in key/value pairs.
#
# key => filename ( e.g. filename.suffix )
# value => absolue filepath ( e.g. /path/to/filename.suffix )
#
# $file_count - The number of files stored with file extension *.suffix.
sub find_files{
    
	# The input file directory.
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	# The file extension suffix.
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
	if(-d $infile_dir){ # Check if $infile_dir is a directory.
		my %files = ();
		my $file_counter = 0;
		opendir(DIR, $infile_dir) || die "Error in opening dir $infile_dir\n";
		while( my $file_name = readdir(DIR)){
			my $infile_name = join('/', $infile_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
			warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
			$files{$file_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
			$file_counter++ if ($file_name =~ m/\.$suffix$/);
		}
		closedir(DIR);
		
		return (\%files, $file_counter);
	}else{
		die "Error $infile_dir does not exist!\n";
	}
}

# sub maintenance_start{	
# 
# 	my $website_name = shift;
# 	die "Error lost input file directory" unless defined $website_name;
# 	my $webapp_dir = shift;
# 	die "Error lost input file directory" unless defined $webapp_dir;
# 
# 	warn "Starting maintenance for $website_name....\n";
# 	my $outfile = join('/', $webapp_dir, "tmp", "maintenance.html");
# 	open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
# 	print OUTFILE <<"EOF";
#     
# 	<html>
# 		<head>
# 		<title>$webapp_name is temporarily unavailable</title>
# 		<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
# 		<style>
# 			body { text-align: center; padding: 150px; }
# 			h1 { font-size: 50px; }
# 			body { font: 20px Helvetica, sans-serif; color: #333; }
# 			article { display: block; text-align: center; width: 650px; margin: 0 auto; }
# 		</style>
# 		</head>
# 		<body>
# 			<article>
# 				<h1>$webapp_name is temporarily unavailable.</h1>
# 				<p>We are currently performing a bulk upload of LIMS samples and will be back online shortly.</p>
# 			</article>
# 		</body>
# 	</html>
#     
# EOF
# 	close(OUTFILE) or die "Couldn't close file $outfile"; 
# }

# sub maintenance_stop{
# 
# 	my $website_name = shift;
# 	die "Error lost the website name" unless defined $website_name;
# 
# 	warn "Stopping maintenance for $website_name....\n";
# 
# 	my $outfile = join('/', $webapp_dir, "tmp", "maintenance.html");
# 	unlink $outfile or warn "Could not unlink $outfile: $!";
# }