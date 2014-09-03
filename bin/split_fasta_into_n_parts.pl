#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
use POSIX;
# perl split_fasta_into_n_parts.pl -i /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_db_2014-05-30/ncbi_nr_db-2014-05-30.fasta -s 3GB -o ~/workspace/testing

my ($fasta_infile, $file_part_size, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      's=s'    => \$file_part_size,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $fasta_infile
      and defined $output_dir
);

$file_part_size = '5GB' unless defined $file_part_size;

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -s file_part_size -o output_dir
    
Description - 
    
OPTIONS:

      -i fasta_infile - 
    
      -s file_part_size - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}


my $part_size = 0;
if($file_part_size =~ m/(-*\d+)(B|MB|GB)/){
      my ($num_size, $byte_size) = ($1, $2);
      die "Error: $file_part_size has to be a positive number greater than 0" unless($num_size > 0);
      if($byte_size eq "B"){
	    $part_size = ($num_size * 1024 );
      }elsif($byte_size eq "MB"){
	    $part_size = ($num_size * (1024 ** 2));
      }elsif($byte_size eq "GB"){
	    $part_size = ($num_size * (1024 ** 3));
      }
}else{
      die "Error: file_part_size must be one of XB or XMB or XGB where X is a positive number greater than 0";
}
my $filesize = (-s $fasta_infile);

my $total_file_count = ceil($filesize / $part_size);

# my $fasta_filename = fileparse($fasta_infile);
my($fasta_filename, $fasta_dir, $fasta_suffix) = fileparse($fasta_infile, qr/\.[^.]*/);
warn "Splitting $fasta_filename into $total_file_count files of size $file_part_size....\n\n";

my $part = 1;
my $size = 0;
my $fasta_outfile = join('/', $output_dir, join("", $fasta_filename, ".part", $part, $fasta_suffix));
warn join("", $fasta_filename, ".part", $part, $fasta_suffix) . "\n";
open(OUTFILE, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!"; 
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $seq_desc = $seq_entry->desc;

      my $fasta_unique_id = $sequence;
      $fasta_unique_id = join(" ", $seq_id, $seq_desc) if(defined($seq_entry->desc));
# 	    warn $seq_id . "\n";
# 	    warn $sequence . "\n";
      my $fasta_seq = join("\n", ">$fasta_unique_id", $sequence);
      print OUTFILE $fasta_seq . "\n";
      $size += length $fasta_seq;
      if($size >= $part_size){
	    warn "File size $size" . "\n";
	    close(OUTFILE) or die "Couldn't close file $fasta_outfile";
	    $part++;
	    my $fasta_outfile = join('/', $output_dir, join("", $fasta_filename, ".part", $part, $fasta_suffix));
	    warn join("", $fasta_filename, ".part", $part, $fasta_suffix) . "\n";
	    open(OUTFILE, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!"; 
	    $size = 0;
      }
}

# Clean out the sequence I/O object.
$seqio = ();

