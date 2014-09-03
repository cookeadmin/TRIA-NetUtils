#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;
use DBD::Pg;
use Switch;
use open qw/:std :utf8/;

# perl postgres_dbquery_writer.pl -i /TRIA-NetUtils/reference_lists/cooke-db-list -o ~/workspace/databases
my ($db_file, $output_dir);
GetOptions(
	'i=s'    => \$db_file,
	'o=s'    => \$output_dir
);

usage() unless (
      defined $output_dir
);

$db_file = "$ENV{\"TRIANetUtils\"}/reference_lists/cooke-db-list" unless defined $db_file;

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

my ($db_host, $db_super_user, $db_password);
$db_host = "localhost";
$db_super_user = "postgres";
$db_password = "postgresql*ygg";
warn $db_file . "\n";

my $db_name_list = get_db_names($db_file);
foreach my $db_name (@{$db_name_list}){

      warn $db_name . "\n";

      my $dbh = DBI->connect(
	    "dbi:Pg:dbname=$db_name;host=$db_host", 
	    $db_super_user, 
	    $db_password
      ) or die "Unable to connect: $DBI::errstr\n";

      # Grabbing public databases
      my $sth = $dbh->prepare("
	    SELECT 'postgresql' AS dbms, t.table_catalog, t.table_schema,t.table_name, 
	    c.column_name, c.ordinal_position, c.data_type, c.character_maximum_length, 
	    n.constraint_type, k2.table_schema, k2.table_name, k2.column_name 
	    FROM information_schema.tables t 
	    NATURAL LEFT JOIN information_schema.columns c 
	    LEFT JOIN(information_schema.key_column_usage k 
	    NATURAL JOIN information_schema.table_constraints n 
	    NATURAL LEFT JOIN information_schema.referential_constraints r) 
	    ON c.table_catalog=k.table_catalog 
	    AND c.table_schema=k.table_schema 
	    AND c.table_name=k.table_name 
	    AND c.column_name=k.column_name 
	    LEFT JOIN information_schema.key_column_usage k2 
	    ON k.position_in_unique_constraint=k2.ordinal_position 
	    AND r.unique_constraint_catalog=k2.constraint_catalog 
	    AND r.unique_constraint_schema=k2.constraint_schema 
	    AND r.unique_constraint_name=k2.constraint_name 
	    WHERE t.table_type='BASE TABLE' 
	    AND t.table_schema 
	    NOT IN('information_schema','pg_catalog')
      ") or die "Cannot prepare statement: $DBI::errstr\n";

      $sth->execute;


      my @table_entries;
      while(my @row = $sth->fetchrow()){
	    #warn join(",", @row) . "\n";
	    for(my $i = 0; $i <= $#row; $i++){
		  unless(defined($row[$i])){
			$row[$i] = "";
		  }
		  if(defined($row[$i])){
			if(length($row[$i]) == 0){
			      $row[$i] = "";
			}
			if($row[$i] eq ""){
			      $row[$i] = "";
			}
			$row[$i] =~ s/\n|\r//g;
		  }
	    }
	    push(@table_entries, join("\t", @row));
      }
      $sth->finish; 

      $db_name =~ s/_/-/g;
      my $db_outfile = join('/', $output_dir, join("-", $db_name, "ER-diagram-info.txt"));

      open OUTFILE, ">$db_outfile" or die "Error opening $db_outfile for writting: $!";
      print OUTFILE join("\n", @table_entries) . "\n";
      close OUTFILE or die "Error closing $db_outfile: $!";

      my $rc = $dbh->disconnect() or die "Unable to disconnect: $DBI::errstr\n";
}


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
