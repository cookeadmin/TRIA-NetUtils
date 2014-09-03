#!/usr/bin/perl
use warnings;
use strict;

use LWP::Simple;
use Getopt::Long;

# perl get_gb_flat_files.pl -i /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_organisms_db_2014-06-03/ncbi_nr_db_2014-05-30_seq_descs.txt -d protein -o ~/workspace/adriana/Phloem/phloem-microarray-annotations/ncbi_nr_flat_files   
my ($id_list_file, $database, $output_dir);
GetOptions(
      'i=s'    => \$id_list_file,
      'd=s'    => \$database,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $id_list_file
      and defined $output_dir
);

sub usage {

die <<"USAGE";

Usage: $0 -i id_list_file -d database -o output_dir

Description - 

OPTIONS:

      -i id_list_file -

      -d database - 

      -o output_dir -

USAGE
}

# Download gene records linked to a set of proteins corresponding to a list
# of GI numbers.
$database = 'protein' unless defined $database;

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

#input UIDs in $db (protein GIs)
open(INFILE, "<$id_list_file") or die "Couldn't open file $id_list_file for reading, $!";
my %unique_id_list = ();
my $unique_id_counter = 0;
while(<INFILE>){
      chomp $_;
      my $seq_desc = $_;
      if ($seq_desc =~ /gi\|(\d+)\|.+/){
	    my $unique_id = $1;
	    $unique_id_list{$unique_id} = $seq_desc;
# 	    warn $unique_id . "\n";
	    $unique_id_counter++;
      }
}
close(INFILE) or die "Couldn't close file $id_list_file";

my $flat_outfile = join('/', $output_dir, "ncbi_nr_db_flat_file.txt");
unless(-s $flat_outfile){
      warn "There are $unique_id_counter ids to fetch from eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi.\n";
      my @bulk_flat_file = ();
      foreach my $unique_id (sort {$a <=> $b} keys %unique_id_list) {
	    warn "$unique_id\n"; 
	    #assemble the efetch URL
	    my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$database&id=$unique_id&rettype=gb&retmode=text";

	    #post the efetch URL
	    my $flat_file = get($url);
	    my $unique_flat_file = join("\n", join(" ", "DESCRIPTION", $unique_id_list{$unique_id}), $flat_file);

	    push(@bulk_flat_file, $unique_flat_file);
      }

      open(OUTFILE, ">$flat_outfile") or die "Couldn't open file $flat_outfile for writting, $!";                          
      foreach my $flat_file_entry (@bulk_flat_file){
	    print OUTFILE $flat_file_entry . "\n";
      }
      close(OUTFILE) or die "Couldn't close file $flat_outfile";
}


my @flat_file_notes = ();
open(INFILE, "<$flat_outfile") or die "Couldn't open file $flat_outfile for reading, $!";
while(<INFILE>){
      chomp $_;
      push(@flat_file_notes, $_);
}
close(INFILE) or die "Couldn't close file $flat_outfile";

my $flat_notes = join("\n", @flat_file_notes) . "\n";

# Split on "//"
my @split_flat_file = split(/\/\//, $flat_notes);
my %flat_notes = ();
foreach my $flat_file_entry (@split_flat_file){
      my @split_flat_file_entry = split("\n", $flat_file_entry);
      my ($unique_id, $notes) = "";
      my $is_notes_tag = "false";
      my @special_notes = ();
      foreach my $flat_file_line (@split_flat_file_entry){

	    if($flat_file_line =~ m/DESCRIPTION\s(.+)/){
		  $unique_id = $1;
	    }
	    
	    if($flat_file_line =~ m/^\s+\/note="(.+)"$/){
		  $notes = $1;
		  push(@{$flat_notes{$unique_id}}, $notes);
		  next;
	    }

	    if($flat_file_line =~ m/^\s+\/note="(.+)$/){

		  $notes = $1;
		  push(@special_notes, $notes);
		  $is_notes_tag = "true";
		  next;

	    }

	    if($is_notes_tag eq "true" and $flat_file_line =~ m/^\s+(.+)$/ 
		  and $flat_file_line !~ m/^\s+.+"$/ and $flat_file_line !~ m/^\s+\/\w+=.+$/){

		  $notes = $1;
		  push(@special_notes, $notes);
		  next;

	    }
	    if($is_notes_tag eq "true" and $flat_file_line =~ m/^\s+(.+)"$/ 
		  and $flat_file_line !~ m/^\s+\/\w+=.+$/){
		  
		  $notes = $1;
		  push(@special_notes, $notes);
		  $is_notes_tag = "false";
		  next;

	    }
	    
# 	    Try to incorporate the case where the /notes tag wraps around multiple lines in the genbank flat file. done!!!
	    
      }

      if(scalar(@special_notes) > 0){
	    push(@{$flat_notes{$unique_id}}, join(" ", @special_notes));
      }
}

my $flat_notes_outfile = join('/', $output_dir, "ncbi_nr_db_notes.txt");
open(OUTFILE, ">$flat_notes_outfile") or die "Couldn't open file $flat_notes_outfile for writting, $!";
print OUTFILE join("\t", "unique_id", "description") . "\n";
foreach my $unique_id (sort {$a cmp $b} keys %flat_notes){
      if(scalar(@{$flat_notes{$unique_id}}) > 0){
	    my @flat_notes = ();
	    push(@flat_notes, @{$flat_notes{$unique_id}});
	    my %seen = ();
	    my @unique_flat_notes = grep { ! $seen{ $_ }++ } @flat_notes;
	    my @unique_flat_notes_sorted = sort {$a cmp $b} @unique_flat_notes;
	    print OUTFILE join("\t", $unique_id, join("\001", @unique_flat_notes_sorted)) . "\n";
      }
}
close(OUTFILE) or die "Couldn't close file $flat_notes_outfile";