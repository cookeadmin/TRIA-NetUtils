#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;
use DBD::Pg;
use Switch;
use open qw/:std :utf8/;

# perl database_schemas.pl -i /TRIA-NetUtils/reference_lists/cooke-db-list -o ~/workspace/databases
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

	# Create output directory if it doesn't already exist.
	my $db_dir = join('/', $output_dir, $db_name);
	unless(-d $db_dir){
		mkdir($db_dir, 0777) or die "Can't make directory: $!";
	}

	#warn $db_name . "\n";

	my $dbh = DBI->connect(
			"dbi:Pg:dbname=$db_name;host=$db_host", 
			$db_super_user, 
			$db_password
	) or die "Unable to connect: $DBI::errstr\n";

	# Grabbing public databases
	my $sth = $dbh->prepare("
	SELECT table_name
	FROM information_schema.tables
	WHERE table_schema='public'
	ORDER BY table_name") or die "Cannot prepare statement: $DBI::errstr\n";

	$sth->execute;


	my @table_names;
	while(my @row = $sth->fetchrow()){
		#warn join(",", @row) . "\n";
		
	    for(my $i = 0; $i <= $#row; $i++){
		$row[$i] =~ s/\r\n/; /g;
		$row[$i] =~ s/; +$//;
		$row[$i] =~ s/,\s*/, /g;
		$row[$i] =~ s/n\/a/N\/A/g;
		$row[$i] =~ s/\t/ /g;
		$row[$i] =~ s/\s{2,}/ /g;
		$row[$i] =~ s/^\s+|\s+$//g;
		  
	    }
	    push(@table_names, @row);
	}
	$sth->finish; 
	
	my $table_dir = join('/', $db_dir, "tables");
	unless(-d $table_dir){
		mkdir($table_dir, 0777) or die "Can't make directory: $!";
	}

	warn "\nDownloading the following tables from $db_name.....\n\n";
	foreach my $table (@table_names){
		warn "$table\n";
		$sth = $dbh->prepare("
		SELECT column_name 
		FROM information_schema.columns 
		WHERE table_name = '$table'") or die "Cannot prepare statement: $DBI::errstr\n";
		$sth->execute;

		my @attribute_names;
		my @row = ();
		while(@row = $sth->fetchrow()){
			#warn join(" ", $table, @row) . "\n";
			for(my $i = 0; $i <= $#row; $i++){
		$row[$i] =~ s/\r\n/; /g;
		$row[$i] =~ s/; +$//;
		$row[$i] =~ s/,\s*/, /g;
				$row[$i] =~ s/n\/a/N\/A/g;
				$row[$i] =~ s/\t/ /g;
				$row[$i] =~ s/\s{2,}/ /g;
				$row[$i] =~ s/^\s+|\s+$//g;
			}
			push(@attribute_names, @row);
		}
		$sth->finish; 
		#warn join(" ", $table, join(", ", @attribute_names)) . "\n";

		my $attribute_list = join(",", @attribute_names);
		$sth = $dbh->prepare("
		SELECT * 
		FROM $table") or die "Cannot prepare statement: $DBI::errstr\n";

		$sth->execute;

		my @table_entries;
		while(my @row = $sth->fetchrow()){
			for(my $i = 0; $i <= $#row; $i++){
			      unless(defined($row[$i])){
				    $row[$i] = "N/A";
			      }
			      
			      if(defined($row[$i])){
				    if(length($row[$i]) == 0){
				    
					  $row[$i] = "N/A";
				    }
				    if($row[$i] eq ""){
					  $row[$i] = "N/A";
				    }
		$row[$i] =~ s/\r\n/; /g;
		$row[$i] =~ s/; +$//;
		$row[$i] =~ s/,\s*/, /g;
					$row[$i] =~ s/n\/a/N\/A/g;
					$row[$i] =~ s/\t/ /g;
					$row[$i] =~ s/\s{2,}/ /g;
					$row[$i] =~ s/^\s+|\s+$//g;
			      }
			}
			#	warn join("\t", @row) . "\n";
			push(@table_entries, join("\t", @row));


		}
		$sth->finish; 


		my $db_outfile = join('/', $table_dir, $table . ".txt");

		open OUTFILE, ">$db_outfile" or die "Error opening $db_outfile for writting: $!";
		print OUTFILE join("\t", @attribute_names) . "\n";
		print OUTFILE join("\n", @table_entries) . "\n";
		close OUTFILE or die "Error closing $db_outfile: $!";

	}

	# Grabbing sequence table names.
	$sth = $dbh->prepare("
	SELECT sequence_name 
	FROM information_schema.sequences 
	WHERE sequence_schema='public' 
	ORDER BY sequence_name") or die "Cannot prepare statement: $DBI::errstr\n";

	$sth->execute;


	my @seq_id_names;
	while(my @row = $sth->fetchrow()){
#		warn join(",", @row) . "\n";
		push(@seq_id_names, @row);
	}
	$sth->finish; 

	my $seq_dir = join('/', $db_dir, "sequence");
	unless(-d $seq_dir){
		mkdir($seq_dir, 0777) or die "Can't make directory: $!";
	}	

	warn "\nDownloading the following sequence id tables from $db_name.....\n\n";
	foreach my $seq_id (@seq_id_names){
		warn "$seq_id\n";
		$sth = $dbh->prepare("
		SELECT *
		FROM $seq_id") or die "Cannot prepare statement: $DBI::errstr\n";

		$sth->execute;
		my @attribute_keys = ("sequence_name", "last_value", "start_value", "increment_by", 
					"max_value", "min_value", "cache_value", "log_cnt", "is_cycled", 
					"is_called");
		my @attribute_values;
		my $row = $sth->fetchrow_hashref;
		foreach my $keys (@attribute_keys){
			my $values = $row->{$keys};
			if($keys =~ m/is_cycled|is_called/){
				switch ($values) {
					case 1	{ $values = "t" }
					case 0	{ $values = "f" }
					else { print "Default: previous case not true" }
				}
			}
			push(@attribute_values, $values);
			
		}

		my $db_outfile = join('/', $seq_dir, $seq_id . ".txt");
		open OUTFILE, ">$db_outfile" or die "Error opening $db_outfile for writting: $!";
		print OUTFILE join("\t", @attribute_keys) . "\n";
		print OUTFILE join("\t", @attribute_values) . "\n";
		close OUTFILE or die "Error closing $db_outfile: $!";
		$sth->finish; 




	}

	# treelimbs_production and triad_production have a lims namespace and needs special attention to retrieve all the tables. 
        if($db_name eq "treelimbs_production" or $db_name eq "triad_production" or $db_name eq "treelimbs_test" or $db_name eq "triad_demo_production" or $db_name eq "tria_lims_production"){
	    # Grabbing lims databases
	    my $sth = $dbh->prepare("
	    SELECT table_name
	    FROM information_schema.tables
	    WHERE table_schema='lims'
	    ORDER BY table_name") or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;


	    my @table_names;
	    while(my @row = $sth->fetchrow()){
		  for(my $i = 0; $i <= $#row; $i++){
		$row[$i] =~ s/\r\n/; /g;
		$row[$i] =~ s/; +$//;
		$row[$i] =~ s/,\s*/, /g;
				$row[$i] =~ s/n\/a/N\/A/g;
				$row[$i] =~ s/\t/ /g;
				$row[$i] =~ s/\s{2,}/ /g;
				$row[$i] =~ s/^\s+|\s+$//g;
		  }
		  #warn join(",", @row) . "\n";
		  push(@table_names, @row);
	    }
	    $sth->finish;

	    warn "\nDownloading the following tables from $db_name.....\n\n";
	    foreach my $table (@table_names){
		  warn "$table\n";
		  $sth = $dbh->prepare("
		  SELECT column_name 
		  FROM information_schema.columns 
		  WHERE table_name = '$table\'") or die "Cannot prepare statement: $DBI::errstr\n";
		  $sth->execute;

		  my @attribute_names;
		  my @row = ();
		  while(@row = $sth->fetchrow()){

			for(my $i = 0; $i <= $#row; $i++){
		$row[$i] =~ s/\r\n/; /g;
		$row[$i] =~ s/; +$//;
		$row[$i] =~ s/,\s*/, /g;
				$row[$i] =~ s/n\/a/N\/A/g;
				$row[$i] =~ s/\t/ /g;
				$row[$i] =~ s/\s{2,}/ /g;
				$row[$i] =~ s/^\s+|\s+$//g;
			}
			#warn join(" ", $table, @row) . "\n";
			push(@attribute_names, @row);
		  }
		  $sth->finish; 
		  #warn join(" ", $table, join(", ", @attribute_names)) . "\n";

		  my $attribute_list = join(",", @attribute_names);
		  $sth = $dbh->prepare("
		  SELECT * 
		  FROM lims.$table") or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;

		  my @table_entries;
		  while(my @row = $sth->fetchrow()){
			for(my $i = 0; $i <= $#row; $i++){
			      unless(defined($row[$i])){
				    $row[$i] = "N/A";
			      }
			      
			      if(defined($row[$i])){
				if(length($row[$i]) == 0){

					$row[$i] = "N/A";
				}
				if($row[$i] eq ""){
					$row[$i] = "N/A";
				}
		$row[$i] =~ s/\r\n/; /g;
		$row[$i] =~ s/; +$//;
		$row[$i] =~ s/,\s*/, /g;
				$row[$i] =~ s/n\/a/N\/A/g;
				$row[$i] =~ s/\t/ /g;
				$row[$i] =~ s/\s{2,}/ /g;
				$row[$i] =~ s/^\s+|\s+$//g;
			      }
			}
	    #	warn join("\t", @row) . "\n";
			push(@table_entries, join("\t", @row));
			      

		  }
		  $sth->finish; 


		  my $db_outfile = join('/', $table_dir, $table . ".txt");

		  open OUTFILE, ">$db_outfile" or die "Error opening $db_outfile for writting: $!";
		  print OUTFILE join("\t", @attribute_names) . "\n";
		  print OUTFILE join("\n", @table_entries) . "\n";
		  close OUTFILE or die "Error closing $db_outfile: $!";
		  
	    }

	    # Grabbing sequence table names.
	    $sth = $dbh->prepare("
	    SELECT sequence_name 
	    FROM information_schema.sequences 
	    WHERE sequence_schema='lims' 
	    ORDER BY sequence_name") or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;


	    my @seq_id_names;
	    while(my @row = $sth->fetchrow()){
      #		warn join(",", @row) . "\n";
		  push(@seq_id_names, @row);
	    }
	    $sth->finish; 

	    my $seq_dir = join('/', $db_dir, "sequence");
	    unless(-d $seq_dir){
		  mkdir($seq_dir, 0777) or die "Can't make directory: $!";
	    }	

	    warn "\nDownloading the following sequence id tables from $db_name.....\n\n";
	    foreach my $seq_id (@seq_id_names){
		  warn "$seq_id\n";
		  $sth = $dbh->prepare("
		  SELECT *
		  FROM lims.$seq_id") or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;
		  my @attribute_keys = ("sequence_name", "last_value", "start_value", "increment_by", 
					  "max_value", "min_value", "cache_value", "log_cnt", "is_cycled", 
					  "is_called");
		  my @attribute_values;
		  my $row = $sth->fetchrow_hashref;
		  foreach my $keys (@attribute_keys){
			      my $values = $row->{$keys};
			      if($keys =~ m/is_cycled|is_called/){
				    switch ($values) {
					  case 1	{ $values = "t" }
					  case 0	{ $values = "f" }
					  else { print "Default: previous case not true" }
				    }
			      }
			      push(@attribute_values, $values);
			      
		  }

		  my $db_outfile = join('/', $seq_dir, $seq_id . ".txt");
		  open OUTFILE, ">$db_outfile" or die "Error opening $db_outfile for writting: $!";
		  print OUTFILE join("\t", @attribute_keys) . "\n";
		  print OUTFILE join("\t", @attribute_values) . "\n";
		  close OUTFILE or die "Error closing $db_outfile: $!";
		  $sth->finish; 




	    }
	}
	
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
