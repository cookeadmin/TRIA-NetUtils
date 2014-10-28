#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;
use DBD::Pg;
use IPC::Open3;
use Spreadsheet::ParseExcel;
use POSIX;
use open qw/:std :utf8/;

#  sudo perl lims_database_bulk_upload.pl -i ~/workspace/bulk_upload_files/reformated_csv_files/M004_fungal_associates.csv -n triad_demo -o ~/workspace/triad_lims_demo_backups
my ($infile, $webapp_name, $backup_dir);
GetOptions(
	'i=s'    => \$infile,
	'n=s'    => \$webapp_name,
	'o=s'    => \$backup_dir
);

usage() unless (
      defined $infile
      and defined $webapp_name
      and defined $backup_dir
);

my ($db_name, $webapp_dir);
if($webapp_name eq 'TRIA-LIMS'){

      $db_name = 'triad_production';
      $webapp_dir = '/var/www/apps/triad';

}elsif($webapp_name eq 'triad-demo'){

      $db_name = 'tria_lims_production';
      $webapp_dir = '/var/www/apps/triad_demo';

}else{
      warn "\nYou entered $webapp_name. Please re-enter the -n webapp_name parameter (must be triad-demo or TRIA-LIMS).\n";
      usage();
}

my $pg_dump	= '/usr/bin/pg_dump';

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i infile -n webapp_name -o backup_dir
    
Description - 
    
OPTIONS:
      -i infile - 
      -n webapp_name - 
      -o backup_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $backup_dir){
	mkdir($backup_dir, 0777) or die "Can't make directory: $!";
}

my ($db_host, $db_super_user, $db_password);
$db_host = "localhost";
$db_super_user = "postgres";
$db_password = "postgresql*ygg";

warn "Preparing some_file for bulk upload into $db_name......" . "\n\n";

warn "Dumping database $db_name to $backup_dir......\n";

my $date_stamp = strftime("%Y-%m-%d_%H-%M-%S", localtime);
# my $date_stamp = strftime("%Y-%m-%d", localtime);

$ENV{'PGPASSWORD'} = $db_password;

my $sqlfilename = join("_", $db_name, $date_stamp) . ".psql";
my $sqlfile = join('/', $backup_dir, $sqlfilename);
warn $sqlfile . "\n\n";

system($pg_dump,
      '--clean',
      '-U', $db_super_user,
      $db_name,
      '-w', 
      '-f', $sqlfile
) == 0 or die "Error executing pg_dump: $?";

maintenance_start($webapp_name, $webapp_dir);

warn "Connecting to database $db_name\n";

# initialize database handlers.
my($dbh, $sth);

$dbh = DBI->connect(
      "dbi:Pg:dbname=$db_name;host=$db_host", 
      $db_super_user, 
      $db_password
) or die "Unable to connect: $DBI::errstr\n";


open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
my $i = 0;
my (%column_index, %column_values);
my $row_counter = 0;
while(<INFILE>){
      chomp $_;
      my @split_entry = split(/\t/, $_);
      if($i eq 0){
	    my $j = 0;
	    foreach my $column_name (@split_entry){
		  $column_index{$j} = $column_name;
		  $j++;
	    }
      }
      else{
	    my $j = 0;
	    foreach my $column_value (@split_entry){
		  push(@{$column_values{$column_index{$j}}}, $column_value);
		  $j++;
	    }
	    $row_counter++;
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $infile";

for(my $i = 0; $i < $row_counter; $i++){

# 	    id	landscape_name	stand_name	code	field_code	experiment_id	comment
      my $unique_sample_id = @{$column_values{"Unique Sample ID"}}[$i] if(defined(@{$column_values{"Unique Sample ID"}}[$i])) or die "Error: Unique Sample ID is blank on line $i of file $infile";

      # Check if the unique sample id matches the sample code critera.
      # ex. M004-13-01-01ST-UL01DA01G01
      ($unique_sample_id =~ m/\w\d+-\d+-\d+-\d+\w*-\w+\d*/) or die "Error: The unique sample id: $unique_sample_id is not the proper format for uploading into the database.\nPlease use the following format; M000-00-00-tree_code-organism_code";
		  
      # Split the unique sample id into experiment, landscape, stand, tree, and organism codes.
      my @split_unique_sample_id = split(/-/, $unique_sample_id);

      my ($experiment_sample_code, $landscape_sample_code, $stand_sample_code, $tree_sample_code, $organism_sample_code) = @split_unique_sample_id;
      warn "Processing Unique Sample ID: " . join("-", $experiment_sample_code, $landscape_sample_code, $stand_sample_code, $tree_sample_code, $organism_sample_code) . "\n";

      # Grab the experiment id number so that we can link tables back to the experiment code.
      $sth = $dbh->prepare("
	    SELECT id
	    FROM lims.experiments
	    WHERE code='$experiment_sample_code'
      ") or die "Cannot prepare statement: $DBI::errstr\n";

      $sth->execute;

      my $experiment_database_id = $sth->fetchrow() or die "Error: There is no experiment id associated with $experiment_sample_code in database $db_name.\nPlease assign yourself a new experiment id using lims new experiment.";

      $sth->finish;
      
      # Checking experiment id.
      warn "experiment_database_id <=> $experiment_database_id\n";

      my ($landscape_id,$landscape_desc_field_code,$stand_id,$stand_desc_name,$landscape_stand_comments,$landscape_name,$stand_name,$landscape_stand_code,$stand_field_code,$experiment_id,$comment_landscape_stand);

      # Grabbing landscape and stand information for populating the landscape_stands table.
      # Landscape ID	Landscape Description - Field Code	Stand ID	Stand Description - Name
      $landscape_id = @{$column_values{"Landscape ID"}}[$i] if(defined(@{$column_values{"Landscape ID"}}[$i])) or die "Error: Landscape ID is blank on line $i of file $infile";

      $landscape_desc_field_code = @{$column_values{"Landscape Description - Field Code"}}[$i] if(defined(@{$column_values{"Landscape Description - Field Code"}}[$i])) or die "Error: Landscape Description - Field Code is blank on line $i of file $infile";

      $stand_id = @{$column_values{"Stand ID"}}[$i] if(defined(@{$column_values{"Stand ID"}}[$i])) or die "Error: Stand ID is blank on line $i of file $infile";
      $stand_id = "0" . $stand_id if(($stand_id >= 1 and $stand_id < 10) and $stand_id !~ m/0\d+/);
      
      $stand_desc_name = @{$column_values{"Stand Description - Name"}}[$i] if(defined(@{$column_values{"Stand Description - Name"}}[$i])) or die "Error: Stand Description - Name is blank on line $i of file $infile";

      $landscape_stand_comments = @{$column_values{"Landscape Stand - Comments"}}[$i] if(defined(@{$column_values{"Landscape Stand - Comments"}}[$i])) or die "Error: Landscape Stand - Comments is blank on line $i of file $infile";

      # Checking the contents from the spreadsheet for landscape and stand data.
      warn join(", ", "Landscape ID", "Landscape Description - Field Code", "Stand ID","Stand Description - Name") . "\n";
      warn join(", ", $landscape_id, $landscape_desc_field_code, $stand_id, $stand_desc_name) . "\n";
      

      # landscape_name	stand_name	code	field_code	experiment_id	comment
      $landscape_name = $landscape_desc_field_code;

      $stand_name = $stand_desc_name;
      
      # Change $infile in ERROR message to file basename when forloop for files are working.
      # Assign $landscape_stand_code by joining the $landscape_sample_code and $stand_sample_code together if file columns match.
      $landscape_stand_code = join("-", $landscape_sample_code, $stand_sample_code) if(($landscape_id eq $landscape_sample_code) and ($stand_id eq $stand_sample_code)) or die "Error: Landscape and Stand from file $infile doesn't match the unique sample id Landscape and Stand codes\nLandscape: $landscape_id <=> $landscape_sample_code\nStand: $stand_id <=> $stand_sample_code";

      # Assign $field code if Stand Description - Field Code is present in the file. Otherwise give it a unique field code by taking the first three letters of $landscape_desc_field_code and the $stand_id number.
      $stand_field_code = @{$column_values{"Stand Description - Field Code"}}[$i] if(defined(@{$column_values{"Stand Description - Field Code"}}[$i]));
      $stand_field_code = join("-", substr($landscape_desc_field_code, 0, 3), sprintf("%d", $stand_id)) unless(defined($stand_field_code));

      $experiment_id = $experiment_database_id;

      $comment_landscape_stand = $landscape_stand_comments;

      warn join(", ", $landscape_name, $stand_name, $landscape_stand_code, $stand_field_code, $experiment_id, $comment_landscape_stand)  . "\n";

      # Grab the landscape_stand id number so that we can link tables back to the landscape_stand entry.
      $sth = $dbh->prepare("
	    SELECT id
	    FROM lims.landscape_stands
	    WHERE landscape_name='$landscape_name'
	    AND stand_name='$stand_name'
	    AND code='$landscape_stand_code'
	    AND field_code='$stand_field_code'
	    AND experiment_id='$experiment_id'
      ") or die "Cannot prepare statement: $DBI::errstr\n";

      $sth->execute;

      # Checking if the select query returned nothing so we give it the value "NULL".
      my $landscape_stand_database_id = $sth->fetchrow();
      if(!defined($landscape_stand_database_id)){
	    $landscape_stand_database_id = "NULL";

      }
      $sth->finish;

      # If the landscape stand id from the database is "NULL" that means that we need to populate the database with a new landscape stand entry.
      if($landscape_stand_database_id =~ m/NULL/){
	    my $landscape_stand_insert_command = "INSERT INTO lims.landscape_stands (landscape_name,stand_name,code,field_code,experiment_id,comment) VALUES ('$landscape_name','$stand_name','$landscape_stand_code','$stand_field_code','$experiment_id','$comment_landscape_stand')";
	    warn $landscape_stand_insert_command . "\n";
	    $sth = $dbh->prepare(
		  $landscape_stand_insert_command
	    ) or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;
	    $sth->finish;

	    # Checking if the insert command worked or not.
	    $sth = $dbh->prepare("
		  SELECT COUNT(*)
		  FROM lims.landscape_stands
		  WHERE landscape_name='$landscape_name'
		  AND stand_name='$stand_name'
		  AND code='$landscape_stand_code'
		  AND field_code='$stand_field_code'
		  AND experiment_id='$experiment_id'
	    ") or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;
	    my $landscape_stand_id_count = $sth->fetchrow();
	    die "Error: $landscape_stand_insert_command\nThe insert command accidentally added another landscape stand id. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($landscape_stand_id_count > 1);
	    die "Error: $landscape_stand_insert_command\nThe insert command did not add the landscape stand entry correctly. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($landscape_stand_id_count eq 0);
	    $sth->finish;

	    # If the insert command worked then grab the landscape_stand id number so that we can link tables back to the landscape_stand entry.
	    if($landscape_stand_id_count eq 1){
		  $sth = $dbh->prepare("
			SELECT id
			FROM lims.landscape_stands
			WHERE landscape_name='$landscape_name'
			AND stand_name='$stand_name'
			AND code='$landscape_stand_code'
			AND field_code='$stand_field_code'
			AND experiment_id='$experiment_id'
		  ") or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;
		  $landscape_stand_database_id = $sth->fetchrow();
		  $sth->finish;
	    }
      }
      
      # Checking if $landscape_stand_database_id is a number or not.
      $landscape_stand_database_id =~ m/\d+/ or die "Error: The landscape stand id; $landscape_stand_database_id from database $db_name is not a number. ";

      warn "landscape_stand_database_id <=> $landscape_stand_database_id\n";

      my ($tree_number,$latitude,$longitude,$tree_height,$altitude,$tree_num_code,$tree_field_code,$attack,$age,$stand_aspect,$stand_slope,$stand_density,$stand_composition,$stand_baited,$success,$description_field,$description_lab,$field_tree_comment,$landscape_stand_id,$alternative_code);

      # Grabbing field tree information for populating the field_trees table.
      # tree_number	latitude	longitude	dbh	altitude	field_code	attack	age	stand_aspect	stand_slope	stand_density	stand_composition	stand_baited	success	description_field	description_lab	comment	landscape_stand_id	alternative_code
      $tree_number = $tree_sample_code;

      $latitude = @{$column_values{"Latitude (Decimal Degrees)"}}[$i] if(defined(@{$column_values{"Latitude (Decimal Degrees)"}}[$i]));
      $latitude = "0.0" unless(defined($latitude));

      $longitude = @{$column_values{"Longitude (Decimal Degrees)"}}[$i] if(defined(@{$column_values{"Longitude (Decimal Degrees)"}}[$i]));
      $longitude = "0.0" unless(defined($longitude));

      $tree_height = @{$column_values{"Tree Description - DBH (cm)"}}[$i] if(defined(@{$column_values{"Tree Description - DBH (cm)"}}[$i]));
      $tree_height = @{$column_values{"Tree height"}}[$i] if(defined(@{$column_values{"Tree height"}}[$i]));
      $tree_height = @{$column_values{"DBH (cm)"}}[$i] if(defined(@{$column_values{"DBH (cm)"}}[$i]));
      $tree_height = "0" unless(defined($tree_height)); 
      
      $altitude = @{$column_values{"Altitude"}}[$i] if(defined(@{$column_values{"Altitude"}}[$i]));
      $altitude = "N/A" unless(defined($altitude)); 

      $tree_num_code = "";
      if($tree_number =~ m/(\d+)(\w+)/){
	    $tree_num_code = join("", sprintf("%d", $1), $2);
      }else{
	    $tree_num_code = $tree_number;
      }

      
      $tree_field_code = join("-", substr($landscape_desc_field_code, 0, 3), sprintf("%d", $stand_id), $tree_num_code) if(($stand_id eq $stand_sample_code) and ($tree_number eq $tree_sample_code)) or die "Error: Stand and Tree Number from file $infile doesn't match the unique sample id Stand and Tree Number codes\nStand: $stand_id <=> $stand_sample_code\nTree Number: $tree_number <=> $tree_sample_code";
      $tree_field_code = @{$column_values{"Tree Field Code"}}[$i] if(defined(@{$column_values{"Tree Field Code"}}[$i]));

      $attack = @{$column_values{"A/U Description"}}[$i] if(defined(@{$column_values{"A/U Description"}}[$i]));
      $attack = @{$column_values{"Affected/Unaffected"}}[$i] if(defined(@{$column_values{"Affected/Unaffected"}}[$i]));
      $attack = 0 unless(defined($attack));
      if($attack =~ m/^A$/ or $attack =~ m/^Affected$/){
	    $attack = 1;
      }elsif($attack =~ m/^U$/ or $attack =~ m/^Unaffected/){
	    $attack = 0;
      }


      $age = @{$column_values{"Notes - Tree Age (years)"}}[$i] if(defined(@{$column_values{"Notes - Tree Age (years)"}}[$i]));
      $age = "unknown" unless(defined($age)); 

      $stand_aspect = @{$column_values{"Stand Description - Aspect"}}[$i] if(defined(@{$column_values{"Stand Description - Aspect"}}[$i]));
      $stand_aspect = "N/A" unless(defined($stand_aspect));

      $stand_slope = @{$column_values{"Stand Description - Slope"}}[$i] if(defined(@{$column_values{"Stand Description - Slope"}}[$i]));
      $stand_slope = "N/A" unless(defined($stand_slope));

      $stand_density = @{$column_values{"Stand Description - Density"}}[$i] if(defined(@{$column_values{"Stand Description - Density"}}[$i]));
      $stand_density = "N/A" unless(defined($stand_density));

      $stand_composition = @{$column_values{"Stand Description - Composition"}}[$i] if(defined(@{$column_values{"Stand Description - Composition"}}[$i]));
      $stand_composition = "N/A" unless(defined($stand_composition));

      $stand_baited = @{$column_values{"Stand Description - Baited?"}}[$i] if(defined(@{$column_values{"Stand Description - Baited?"}}[$i]));
      $stand_baited = "N/A" unless(defined($stand_baited));

      $success = @{$column_values{"Tree Description - Attack"}}[$i] if(defined(@{$column_values{"Tree Description - Attack"}}[$i]));
      $success = @{$column_values{"Tree Description - Success"}}[$i] if(defined(@{$column_values{"Tree Description - Success"}}[$i]));
      $success = "N/A" unless(defined($success));

      $description_field = @{$column_values{"Tree Description - Other"}}[$i] if(defined(@{$column_values{"Tree Description - Other"}}[$i]));
      $description_field = @{$column_values{"Tree Description - Field"}}[$i] if(defined(@{$column_values{"Tree Description - Field"}}[$i]));
      $description_field = "N/A" unless(defined($description_field)); 

      $description_lab = @{$column_values{"Tree Description - Lab"}}[$i] if(defined(@{$column_values{"Tree Description - Lab"}}[$i]));
      $description_lab = "N/A" unless(defined($description_lab));

      my @field_tree_comments = ();
      push(@field_tree_comments, @{$column_values{"Field Tree - Comments"}}[$i]) if(defined(@{$column_values{"Field Tree - Comments"}}[$i]) and @{$column_values{"Field Tree - Comments"}}[$i] ne "");
      push(@field_tree_comments, join(" ", "Notes - Other lab notes:", @{$column_values{"Notes - Other lab notes"}}[$i])) if(defined(@{$column_values{"Notes - Other lab notes"}}[$i]) and @{$column_values{"Notes - Other lab notes"}}[$i] ne "");
      push(@field_tree_comments, join(" ", "Notes - Other:", @{$column_values{"Notes - Other"}}[$i])) if(defined(@{$column_values{"Notes - Other"}}[$i]) and @{$column_values{"Notes - Other"}}[$i] ne "");
      push(@field_tree_comments, join(" ", "Attack height:", @{$column_values{"Attack height"}}[$i])) if(defined(@{$column_values{"Attack height"}}[$i]) and @{$column_values{"Attack height"}}[$i] ne "");
      push(@field_tree_comments, join(" ", "Section:", @{$column_values{"Section"}}[$i])) if(defined(@{$column_values{"Section"}}[$i]) and @{$column_values{"Section"}}[$i] ne "");
      push(@field_tree_comments, join(" ", "Section Description:", @{$column_values{"Section Description"}}[$i])) if(defined(@{$column_values{"Section Description"}}[$i]) and @{$column_values{"Section Description"}}[$i] ne "");
      
      if(scalar(@field_tree_comments) > 0){
	    $field_tree_comment = join("; ", @field_tree_comments);
      }else{
	    $field_tree_comment = "N/A" unless(defined($field_tree_comment)); 
      }

      $landscape_stand_id = $landscape_stand_database_id;

      $alternative_code = join(" ", join("", uc(substr($landscape_name, 0, 1)), uc(substr($stand_name, 0, 1))), join("-", sprintf("%d", $stand_id), $tree_num_code));

      warn join("\t", $tree_number,$latitude,$longitude,$tree_height,$altitude,$tree_field_code,$attack,$age,$stand_aspect,$stand_slope,$stand_density,$stand_composition,$stand_baited,$success,$description_field,$description_lab,$field_tree_comment,$landscape_stand_id,$alternative_code) . "\n";

      # Grab the field_trees id number so that we can link tables back to the field_trees entry.
      #tree_number	latitude	longitude	dbh	altitude	field_code	attack	age	stand_aspect	stand_slope	stand_density	stand_composition	stand_baited	success	description_field	description_lab	comment	landscape_stand_id	alternative_code
      $sth = $dbh->prepare("
	    SELECT id
	    FROM lims.field_trees
	    WHERE tree_number='$tree_number'
	    AND latitude='$latitude'
	    AND longitude='$longitude'
	    AND dbh='$tree_height'
	    AND altitude='$altitude'
	    AND field_code='$tree_field_code'
	    AND attack='$attack'
	    AND age='$age'
	    AND stand_aspect='$stand_aspect'
	    AND stand_slope='$stand_slope'
	    AND stand_density='$stand_density'
	    AND stand_composition='$stand_composition'
	    AND stand_baited='$stand_baited'
	    AND success='$success'
	    AND description_field='$description_field'
	    AND description_lab='$description_lab'
	    AND landscape_stand_id='$landscape_stand_id'
	    AND alternative_code='$alternative_code'
      ") or die "Cannot prepare statement: $DBI::errstr\n";

      $sth->execute;

      # Checking if the select query returned nothing so we give it the value "NULL".
      my $field_trees_database_id = $sth->fetchrow();
      if(!defined($field_trees_database_id)){
	    $field_trees_database_id = "NULL";

      }
      $sth->finish;
# 
# 	    # If the field trees id from the database is "NULL" that means that we need to populate the database with a new field trees entry.
      if($field_trees_database_id =~ m/NULL/){
	    my $field_trees_insert_command = "INSERT INTO lims.field_trees (tree_number,latitude,longitude,dbh,altitude,field_code,attack,age,stand_aspect,stand_slope,stand_density,stand_composition,stand_baited,success,description_field,description_lab,comment,landscape_stand_id,alternative_code) VALUES ('$tree_number','$latitude','$longitude','$tree_height','$altitude','$tree_field_code','$attack','$age','$stand_aspect','$stand_slope','$stand_density','$stand_composition','$stand_baited','$success','$description_field','$description_lab','$field_tree_comment','$landscape_stand_id','$alternative_code')";
	    warn $field_trees_insert_command . "\n";
	    $sth = $dbh->prepare(
		  $field_trees_insert_command
	    ) or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;
	    $sth->finish;

	    # Checking if the insert command worked or not.
	    $sth = $dbh->prepare("
		  SELECT COUNT(*)
		  FROM lims.field_trees
			WHERE tree_number='$tree_number'
			AND latitude='$latitude'
			AND longitude='$longitude'
			AND dbh='$tree_height'
			AND altitude='$altitude'
			AND field_code='$tree_field_code'
			AND attack='$attack'
			AND age='$age'
			AND stand_aspect='$stand_aspect'
			AND stand_slope='$stand_slope'
			AND stand_density='$stand_density'
			AND stand_composition='$stand_composition'
			AND stand_baited='$stand_baited'
			AND success='$success'
			AND description_field='$description_field'
			AND description_lab='$description_lab'
			AND landscape_stand_id='$landscape_stand_id'
			AND alternative_code='$alternative_code'
	    ") or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;

	    my $field_trees_id_count = $sth->fetchrow();
	    die "Error: $field_trees_database_id\nThe insert command accidentally added another landscape stand id. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($field_trees_id_count > 1);
	    die "Error: $field_trees_database_id\nThe insert command did not add the landscape stand entry correctly. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($field_trees_id_count eq 0);
	    $sth->finish;

	    # If the insert command worked then grab the landscape_stand id number so that we can link tables back to the landscape_stand entry.
	    if($field_trees_id_count eq 1){
			$sth = $dbh->prepare("
			SELECT id
			FROM lims.field_trees
			WHERE tree_number='$tree_number'
			AND latitude='$latitude'
			AND longitude='$longitude'
			AND dbh='$tree_height'
			AND altitude='$altitude'
			AND field_code='$tree_field_code'
			AND attack='$attack'
			AND age='$age'
			AND stand_aspect='$stand_aspect'
			AND stand_slope='$stand_slope'
			AND stand_density='$stand_density'
			AND stand_composition='$stand_composition'
			AND stand_baited='$stand_baited'
			AND success='$success'
			AND description_field='$description_field'
			AND description_lab='$description_lab'
			AND landscape_stand_id='$landscape_stand_id'
			AND alternative_code='$alternative_code'
		  ") or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;
		  $field_trees_database_id = $sth->fetchrow();
		  $sth->finish;
	    }
      }

      # Grabbing sample condition information for populating the sample conditions table.
      # id	code	experiment_factor_detail_id	rank
      my $condition_code = "";
      if($organism_sample_code =~ m/^(U[A-Z]+)\d+D[A-Z]+\d+G\d+/){
	    $condition_code = $1;
      }elsif($organism_sample_code =~ m/^(D[A-Z]+)\d+G\d+/){
	    $condition_code = $1; 
      }elsif($organism_sample_code =~ m/^(F)\d*/){
	    $condition_code = $1; 
      }elsif($organism_sample_code =~ m/^(M)\d*/){
	    $condition_code = $1; 
      }elsif($organism_sample_code =~ m/^(G)\d*/){
	    $condition_code = $1; 
      }

      my $sample_condition_detail_code = join("-", $experiment_sample_code, $landscape_sample_code, $stand_sample_code, $tree_sample_code, $condition_code);

      my %field_collection_codes = (
	    "DL" => 501,
	    "UC" => 502,
	    "UL" => 503,
	    "UM" => 504,
	    "DA" => 505,
	    "UR" => 506,
	    "UCG" => 507,
	    "UCL" => 508,
	    "UCM" => 509,
	    "UCMG" => 510,
	    "ULM" => 511,
	    "UMG" => 512,
	    "M" => 499,
	    "G" => 500,
      );

      my $biomaterial_species = @{$column_values{"Biomaterial - Species"}}[$i] if(defined(@{$column_values{"Biomaterial - Species"}}[$i])) or die "Error: Biomaterial - Species on line $i of file $infile"; #M038

      my %biomaterial_codes = (
	    "Pinus contorta" => 495,#M004
	    "Pinus banksiana" => 496,#M038
	    "Pinus contorta var. latifolia" => 497,
	    "Pinus contorta var. contorta" => 498,
      );
      warn "$condition_code, $field_collection_codes{$condition_code}";
      die "Error: biomaterial_species blank $biomaterial_species, $biomaterial_codes{$biomaterial_species}" if(!defined($biomaterial_codes{$biomaterial_species}));
      die "Error: condition_code blank $condition_code, $field_collection_codes{$condition_code}" if(!defined($field_collection_codes{$condition_code}));

      # code	experiment_factor_detail_id	rank
      my @sample_conditions = ();
#       push(@sample_conditions, join("\t", $sample_condition_detail_code, 495, 1)); #M004
      push(@sample_conditions, join("\t", $sample_condition_detail_code, $biomaterial_codes{$biomaterial_species}, 1));
      push(@sample_conditions, join("\t", $sample_condition_detail_code, $field_collection_codes{$condition_code}, 2));

      foreach my $sample_condition (@sample_conditions){
	    my @split_sample_condition = split(/\t/, $sample_condition);
	    my ($sample_condition_code, $experiment_factor_detail_id, $sample_condition_rank) = @split_sample_condition;
	    warn join(", ", "code", "experiment_factor_detail_id", "rank") . "\n";
	    warn join(", ", $sample_condition_code, $experiment_factor_detail_id, $sample_condition_rank) . "\n";
	    
	    # Grab the samples number so that we can link tables back to the tissue sample entry.
	    $sth = $dbh->prepare("
		  SELECT id
		  FROM lims.sample_conditions
		  WHERE code='$sample_condition_code'
		  AND experiment_factor_detail_id='$experiment_factor_detail_id'
		  AND rank='$sample_condition_rank'
	    ") or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;

	    # Checking if the select query returned nothing so we give it the value "NULL".
	    my $sample_conditions_database_id = $sth->fetchrow();
	    
	    if(!defined($sample_conditions_database_id)){
		  $sample_conditions_database_id = "NULL";
		  warn $sample_conditions_database_id . "\n";
	    }
	    $sth->finish;

	    # If the landscape stand id from the database is "NULL" that means that we need to populate the database with a new landscape stand entry.
	    if($sample_conditions_database_id =~ m/NULL/){
		  my $sample_conditions_insert_command = "INSERT INTO lims.sample_conditions (code,experiment_factor_detail_id,rank) VALUES ('$sample_condition_code','$experiment_factor_detail_id','$sample_condition_rank')";
		  warn $sample_conditions_insert_command . "\n";
		  $sth = $dbh->prepare(
			$sample_conditions_insert_command
		  ) or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;
		  $sth->finish;

		  # Checking if the insert command worked or not.
		  $sth = $dbh->prepare("
			SELECT COUNT(*)
			FROM lims.sample_conditions
			WHERE code='$sample_condition_code'
			AND experiment_factor_detail_id='$experiment_factor_detail_id'
			AND rank='$sample_condition_rank'
		  ") or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;

		  my $sample_conditions_id_count = $sth->fetchrow();
		  die "Error: $sample_conditions_insert_command\nThe insert command accidentally added another samples id. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($sample_conditions_id_count > 1);
		  die "Error: $sample_conditions_insert_command\nThe insert command did not add the samples entry correctly. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($sample_conditions_id_count eq 0);
		  $sth->finish;

		  # If the insert command worked then grab the landscape_stand id number so that we can link tables back to the landscape_stand entry.
		  if($sample_conditions_id_count eq 1){
			$sth = $dbh->prepare("
			      SELECT id
			      FROM lims.sample_conditions
			      WHERE code='$sample_condition_code'
			      AND experiment_factor_detail_id='$experiment_factor_detail_id'
			      AND rank='$sample_condition_rank'
			") or die "Cannot prepare statement: $DBI::errstr\n";

			$sth->execute;
			$sample_conditions_database_id = $sth->fetchrow();
			$sth->finish;
		  }
	    }
      }
      # id	date	process_date	people	finish	individual	individual_description	storage	comment	process_detail	code	experiment_id	field_tree_id	sample_condition_detail_code	gallery
      # 25	2008-03-06	2008-03-12	Field/Lab	0	23G40	Origin (wood, G40); associated MPB DA08	N/A	Isolation Date: 2008-03-15Tubed Date: 2008-04-10	10 4"-disks with intact bark, phloem and xylem removed from tree using hole saw and transported on ice to the lab. Flame-sterilized wood sample plated onto 1.5% malt extract agar (MEA) at time of processing.  Fungi subcultured onto MEA as they appeared.  Subcultures allowed to grow to maturity then identified based on cultural and microscopic morphology.  Subcultured 	M002-01-01-04A-UM23G40	7	467	M002-01-01-04A-UM	G40
      my ($date_collected,$processing_date,$personnel,$finished,$individual,$individual_description,$storage,$field_sample_comment,$process_detail,$field_sample_code,$field_trees_id,$gallery);

      $date_collected = @{$column_values{"Date Collected"}}[$i] if(defined(@{$column_values{"Date Collected"}}[$i])) or die "Error: Date Collected is blank on line $i of file $infile";

      $processing_date = @{$column_values{"Processing Date"}}[$i] if(defined(@{$column_values{"Processing Date"}}[$i])) or die "Error: Processing Date is blank on line $i of file $infile";

      $personnel = @{$column_values{"Personnel"}}[$i] if(defined(@{$column_values{"Personnel"}}[$i])) or die "Error: Personnel is blank on line $i of file $infile";

      # Is finished field. Assuming that sample entries are not finished until someone changes it.
      $finished = 0; 

      if($organism_sample_code =~ m/^(U[A-Z]+\d+D[A-Z]+\d+G\d+)/){
	    $individual = $1; 
      }elsif($organism_sample_code =~ m/^(D[A-Z]+\d+G\d+)/){
	    $individual = $1; 
      }elsif($organism_sample_code =~ m/^(F\d*)/){
	    $individual = $1; 
      }elsif($organism_sample_code =~ m/^(M\d*)/){
	    $individual = $1; 
      }elsif($organism_sample_code =~ m/^(G\d*)/){
	    $individual = $1; 
      }

      $storage = "CCIS 5-104" unless(defined(@{$column_values{"Storage"}}[$i]));

      $field_sample_code = $unique_sample_id;

      $field_trees_id = $field_trees_database_id;

      my @process_details = ();
      # Biomaterial - Collection	Biomaterial - Culture
      # Gallery	Gallery Description
      if($condition_code =~ m/^U[A-Z]+/){
		push(@process_details, join(" ", "Gallery Description:", @{$column_values{"Gallery Description"}}[$i])) if(defined(@{$column_values{"Gallery Description"}}[$i]) and @{$column_values{"Gallery Description"}}[$i] ne "");

      }elsif($condition_code =~ m/^D[A-Z]+/){
		push(@process_details, join(" ", "Gallery Description:", @{$column_values{"Gallery Description"}}[$i])) if(defined(@{$column_values{"Gallery Description"}}[$i]) and @{$column_values{"Gallery Description"}}[$i] ne "");

      }elsif($condition_code =~ m/^F/){

      }elsif($condition_code =~ m/^M/){

      }elsif($condition_code =~ m/^G/){

      }

      if(scalar(@process_details) > 0){
	    $process_detail = join("; ", @process_details);
      }else{
	    $process_detail = "N/A" unless(defined($process_detail)); 
      }

      # Sample Processing Date	Fungal Establishment Date	Isolation Date	Tubed Date
      my @field_sample_comments = ();
      if($condition_code =~ m/^U[A-Z]+/){
            push(@field_sample_comments, join(" ", "Processing Date:", @{$column_values{"Processing Date"}}[$i])) if(defined(@{$column_values{"Processing Date"}}[$i]) and @{$column_values{"Processing Date"}}[$i] ne "");
	    push(@field_sample_comments, join(" ", "Fungal Establishment Date:", @{$column_values{"Fungal Establishment Date"}}[$i])) if(defined(@{$column_values{"Fungal Establishment Date"}}[$i]) and @{$column_values{"Fungal Establishment Date"}}[$i] ne "");
	    push(@field_sample_comments, join(" ", "Isolation Date:", @{$column_values{"Isolation Date"}}[$i])) if(defined(@{$column_values{"Isolation Date"}}[$i]) and @{$column_values{"Isolation Date"}}[$i] ne "");
	    push(@field_sample_comments, join(" ", "Tubed Date:", @{$column_values{"Tubed Date"}}[$i])) if(defined(@{$column_values{"Tubed Date"}}[$i]) and @{$column_values{"Tubed Date"}}[$i] ne "");
      }elsif($condition_code =~ m/^D[A-Z]+/){
            push(@field_sample_comments, join(" ", "Processing Date:", @{$column_values{"Processing Date"}}[$i])) if(defined(@{$column_values{"Processing Date"}}[$i]) and @{$column_values{"Processing Date"}}[$i] ne "");
	    push(@field_sample_comments, join(" ", "Notes on individual MPB:", @{$column_values{"Comments - Notes on individual MPB"}}[$i])) if(defined(@{$column_values{"Comments - Notes on individual MPB"}}[$i]) and @{$column_values{"Comments - Notes on individual MPB"}}[$i] ne "");
	    push(@field_sample_comments, join(" ", "Used for DNA:", @{$column_values{"Comments - Used for DNA"}}[$i])) if(defined(@{$column_values{"Comments - Used for DNA"}}[$i]) and @{$column_values{"Comments - Used for DNA"}}[$i] ne "");
	    push(@field_sample_comments, join(" ", "Used for fungal sampling:", @{$column_values{"Comments - Used for fungal sampling"}}[$i])) if(defined(@{$column_values{"Comments - Used for fungal sampling"}}[$i]) and @{$column_values{"Comments - Used for fungal sampling"}}[$i] ne "");
	    
      }elsif($condition_code =~ m/^F/){

      }elsif($condition_code =~ m/^M/){

      }elsif($condition_code =~ m/^G/){

      }

      if(scalar(@field_sample_comments) > 0){
	    $field_sample_comment = join("; ", @field_sample_comments);
      }else{
	    $field_sample_comment = "N/A" unless(defined($field_sample_comment)); 
      }

      $gallery = @{$column_values{"Gallery"}}[$i] if(defined(@{$column_values{"Gallery"}}[$i]));
      if(defined($gallery)){
	    if($gallery !~ m/^G\d+/){
		  if(($gallery < 10) and ($gallery !~ m/^0/)){
			$gallery = join("", 0, $gallery);
		  }
		  $gallery = join("", "G", $gallery);
	    }
      }
      $gallery = "N/A" unless(defined($gallery)); 

      # 	    Associated MPB	Substrate Origin (wood, G18)
      my @individual_descriptions = ();
      if($condition_code =~ m/^U[A-Z]+/){
	    push(@individual_descriptions, join(" ", "Substrate:", join("", "(", join(", ", @{$column_values{"Substrate"}}[$i], $gallery), ")"))) if(defined(@{$column_values{"Substrate"}}[$i]) and @{$column_values{"Substrate"}}[$i] ne "");
	    push(@individual_descriptions, join(" ", "Associated MPB:", @{$column_values{"Associated MPB"}}[$i])) if(defined(@{$column_values{"Associated MPB"}}[$i]) and @{$column_values{"Associated MPB"}}[$i] ne "");
      }elsif($condition_code =~ m/^D[A-Z]+/){
	    push(@individual_descriptions, join(" ", "Associated wood sample:", @{$column_values{"Comments - Associated wood sample"}}[$i])) if(defined(@{$column_values{"Comments - Associated wood sample"}}[$i]) and @{$column_values{"Comments - Associated wood sample"}}[$i] ne "");
      }elsif($condition_code =~ m/^F/){

      }elsif($condition_code =~ m/^M/){
	    push(@individual_descriptions, @{$column_values{"Biomaterial - Description"}}[$i]) if(defined(@{$column_values{"Biomaterial - Description"}}[$i]) and @{$column_values{"Biomaterial - Description"}}[$i] ne "");

      }elsif($condition_code =~ m/^G/){
	    push(@individual_descriptions, @{$column_values{"Biomaterial - Description"}}[$i]) if(defined(@{$column_values{"Biomaterial - Description"}}[$i]) and @{$column_values{"Biomaterial - Description"}}[$i] ne "");
      }

      if(scalar(@individual_descriptions) > 0){
	    $individual_description = join("; ", @individual_descriptions);
      }else{
	    $individual_description = "N/A" unless(defined($individual_description)); 
      }

      warn join(",", $date_collected,$processing_date,$personnel,$finished,$individual,$individual_description,$storage,$field_sample_comment,$process_detail,$field_sample_code,$field_trees_id,$sample_condition_detail_code,$gallery) . "\n";

      # Grab the samples number so that we can link tables back to the tissue sample entry.
      $sth = $dbh->prepare("
	    SELECT id
	    FROM lims.field_samples
	    WHERE date='$date_collected'
	    AND process_date='$processing_date'
	    AND people='$personnel'
	    AND finish='$finished'
	    AND individual='$individual'
	    AND individual_description='$individual_description'
	    AND storage='$storage'
	    AND comment='$field_sample_comment'
	    AND process_detail='$process_detail'
	    AND code='$field_sample_code'
	    AND experiment_id='$experiment_id'
	    AND field_tree_id='$field_trees_id'
	    AND sample_condition_detail_code='$sample_condition_detail_code'
	    AND gallery='$gallery'
      ") or die "Cannot prepare statement: $DBI::errstr\n";

      $sth->execute;

      # Checking if the select query returned nothing so we give it the value "NULL".
      my $field_samples_database_id = $sth->fetchrow();
      
      if(!defined($field_samples_database_id)){
	    $field_samples_database_id = "NULL";
	    warn $field_samples_database_id . "\n";
      }
      $sth->finish;

      # If the landscape stand id from the database is "NULL" that means that we need to populate the database with a new landscape stand entry.
      if($field_samples_database_id =~ m/NULL/){
	    my $field_samples_insert_command = "INSERT INTO lims.field_samples (date,process_date,people,finish,individual,individual_description,storage,comment,process_detail,code,experiment_id,field_tree_id,sample_condition_detail_code,gallery) VALUES ('$date_collected','$processing_date','$personnel','$finished','$individual','$individual_description','$storage','$field_sample_comment','$process_detail','$field_sample_code','$experiment_id','$field_trees_id','$sample_condition_detail_code','$gallery')";
	    warn $field_samples_insert_command . "\n";
	    $sth = $dbh->prepare(
		  $field_samples_insert_command
	    ) or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;
	    $sth->finish;

	    # Checking if the insert command worked or not.
	    $sth = $dbh->prepare("
		  SELECT COUNT(*)
		  FROM lims.field_samples
		  WHERE date='$date_collected'
		  AND process_date='$processing_date'
		  AND people='$personnel'
		  AND finish='$finished'
		  AND individual='$individual'
		  AND individual_description='$individual_description'
		  AND storage='$storage'
		  AND comment='$field_sample_comment'
		  AND process_detail='$process_detail'
		  AND code='$field_sample_code'
		  AND experiment_id='$experiment_id'
		  AND field_tree_id='$field_trees_id'
		  AND sample_condition_detail_code='$sample_condition_detail_code'
		  AND gallery='$gallery'
	    ") or die "Cannot prepare statement: $DBI::errstr\n";

	    $sth->execute;

	    my $field_samples_id_count = $sth->fetchrow();
	    die "Error: $field_samples_insert_command\nThe insert command accidentally added another samples id. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($field_samples_id_count > 1);
	    die "Error: $field_samples_insert_command\nThe insert command did not add the samples entry correctly. You will need to backup the database with $sqlfilename and rerun the bulk upload program" if($field_samples_id_count eq 0);
	    $sth->finish;

	    # If the insert command worked then grab the landscape_stand id number so that we can link tables back to the landscape_stand entry.
	    if($field_samples_id_count eq 1){
		  $sth = $dbh->prepare("
			SELECT id
			FROM lims.field_samples
			WHERE date='$date_collected'
			AND process_date='$processing_date'
			AND people='$personnel'
			AND finish='$finished'
			AND individual='$individual'
			AND individual_description='$individual_description'
			AND storage='$storage'
			AND comment='$field_sample_comment'
			AND process_detail='$process_detail'
			AND code='$field_sample_code'
			AND experiment_id='$experiment_id'
			AND field_tree_id='$field_trees_id'
			AND sample_condition_detail_code='$sample_condition_detail_code'
			AND gallery='$gallery'
		  ") or die "Cannot prepare statement: $DBI::errstr\n";

		  $sth->execute;
		  $field_samples_database_id = $sth->fetchrow();
		  $sth->finish;
	    }
      }
}

warn "Disconnecting from database $db_name\n";
my $rc = $dbh->disconnect() or die "Unable to disconnect: $DBI::errstr\n";
maintenance_stop($webapp_name, $webapp_dir);

sub maintenance_start{	

	my $website_name = shift;
	die "Error lost input file directory" unless defined $website_name;
	my $webapp_dir = shift;
	die "Error lost input file directory" unless defined $webapp_dir;

	warn "Starting maintenance for $website_name....\n";
	my $outfile = join('/', $webapp_dir, "tmp", "maintenance.html");
	open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
	print OUTFILE <<"EOF";
    
	<html>
		<head>
		<title>$webapp_name is temporarily unavailable</title>
		<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
		<style>
			body { text-align: center; padding: 150px; }
			h1 { font-size: 50px; }
			body { font: 20px Helvetica, sans-serif; color: #333; }
			article { display: block; text-align: center; width: 650px; margin: 0 auto; }
		</style>
		</head>
		<body>
			<article>
				<h1>$webapp_name is temporarily unavailable.</h1>
				<p>We are currently performing a bulk upload of LIMS samples and will be back online shortly.</p>
			</article>
		</body>
	</html>
    
EOF
	close(OUTFILE) or die "Couldn't close file $outfile"; 
}

sub maintenance_stop{

	my $website_name = shift;
	die "Error lost the website name" unless defined $website_name;

	warn "Stopping maintenance for $website_name....\n";

	my $outfile = join('/', $webapp_dir, "tmp", "maintenance.html");
	unlink $outfile or warn "Could not unlink $outfile: $!";
}