#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;
use DBD::Pg;
use Switch;

my ($db_file, $db_host, $db_super_user, $db_password);
GetOptions(
'i=s'    => \$db_file,
'h=s'    => \$db_host,
's=s'    => \$db_super_user,
'p=s'    => \$db_password,
);
$db_file = "/home/cookeadmin/Desktop/cooke-db-list" unless defined $db_file;

my $db_name_list = get_db_names($db_file);

foreach my $db_name (@{$db_name_list}){
	warn $db_name . "\n";
}

my $db_name = "poplar_production";
$db_host = "localhost";
$db_super_user = "postgres";
$db_password = "postgresql*ygg";

my $dbh = DBI->connect(
		"dbi:Pg:dbname=$db_name;host=$db_host", 
		$db_super_user, 
		$db_password
) or die "Unable to connect: $DBI::errstr\n";

#my $sth = $dbh->prepare("SELECT id, est_id, number, array_row, array_col, row, col, name, contig, mnc_ident, strong_blast_hit, blast_eval, hit_accession, hit_def FROM array_spots") or die "Cannot prepare statement: $DBI::errstr\n";

my $sth = $dbh->prepare("SELECT count(*) FROM users WHERE login = 'muirheadk'") or die "Cannot prepare statement: $DBI::errstr\n";

$sth->execute;

my ($isUser) = $sth->fetchrow_array();

$sth->finish; 

warn "$isUser\n";

if($isUser eq 0){
	my $sth = $dbh->prepare("CREATE ROLE muirheadk LOGIN PASSWORD 'my_password'") or die "Cannot prepare statement: $DBI::errstr\n";
	$sth->execute;
	#my @row = ();
	#while (@row = $sth->fetchrow()){
	#	warn join("\t", @row) . "\n";
	#}
	$sth->finish;
	my $sth = $dbh->prepare("GRANT postgres TO muirheadk") or die "Cannot prepare statement: $DBI::errstr\n";
	$sth->execute;
	#my @row = ();
	#while (@row = $sth->fetchrow()){
	#	warn join("\t", @row) . "\n";
	#}
	$sth->finish;
}

my $rc = $dbh->disconnect() or die "Unable to disconnect: $DBI::errstr\n";


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
