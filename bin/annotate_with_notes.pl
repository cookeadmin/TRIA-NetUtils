#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;

# perl annotate_with_notes.pl -i ~/workspace/adriana/Phloem/phloem-microarray-annotations/blastx_output/ncbi_nr_blastx -f ~/workspace/adriana/Phloem/phloem-microarray-annotations/ncbi_nr_flat_files/ncbi_nr_db_notes.txt -o ~/workspace/adriana/Phloem/phloem-microarray-annotations/blastx_output/ncbi_nr_blastx/annotations_with_notes


my ($input_dir, $flat_notes_file, $output_dir);
my @options = (
    'i=s'    => \$input_dir,
    'f=s'    => \$flat_notes_file,
    'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
    defined $input_dir
    and $flat_notes_file
    and $output_dir
);

sub usage {
    
die << "USAGE";

    Usage: $0 -i input_dir -f flat_notes_file -o output_dir

USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %flat_notes = ();
open(INFILE, "<$flat_notes_file") or die "Couldn't open file $flat_notes_file for reading, $!";
while(<INFILE>){
      chomp $_;
      
      my @split_notes_entry = split(/\t/, $_);
      my ($unique_id, $notes) = ($split_notes_entry[0], $split_notes_entry[1]);
      warn join("\t", $unique_id, $notes) . "\n";
      $flat_notes{$unique_id} = $notes;
}
close(INFILE) or die "Couldn't close file $flat_notes_file";

my $files = find_files($input_dir, "tsv");

foreach my $file (keys %{$files}){
      my $filename = $files->{$file};
      warn $file . "\n";

      my $outfile = join('/', $output_dir, $file . ".notes");
      open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
      open(INFILE, "<$filename") or die "Couldn't open file $filename for reading, $!";
      my $i = 0;
      while(<INFILE>){
	    chomp $_;
	    if($i eq 0){
		  my $blast_header = $_;
		  my $blast_notes_header = join("\t", $blast_header, "notes");
		  print OUTFILE $blast_notes_header . "\n";
	    }else{
		  warn $_ . "\n";
		  my $blast_entry = $_;
		  my @split_blast_entry = split(/\t/, $blast_entry);
		  my $target_id = $split_blast_entry[1];

		  my $notes = "N/A";
		  if(defined($flat_notes{$target_id})){
			$notes = $flat_notes{$target_id};
		  }
		  my $blast_notes_entry = join("\t", $blast_entry, $notes);
		  print OUTFILE $blast_notes_entry . "\n";
	    }
	    $i++;
      }
      close(INFILE) or die "Couldn't close file $filename";
      close(OUTFILE) or die "Couldn't close file $outfile";

}


sub find_files {
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	my $suffix = shift;
	die "Error lost pattern" unless defined $suffix;
    
	opendir(DIR,$dir) or die "Error opening $dir: $!";

	my %files;
	while(my $file = readdir(DIR)){
	    my $filename = join('/', $dir, $file);
	    $files{$file} = $filename if ($file =~ /\.$suffix$/);
	}
	closedir(DIR) or die "Error closing $dir: $!";
	return \%files;
}