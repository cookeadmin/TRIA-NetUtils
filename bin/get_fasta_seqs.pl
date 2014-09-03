#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
# perl get_fasta_seqs.pl -i /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_db_2014-05-30 -l ~/workspace/adriana/organisms-list.txt -o /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_organisms_db

my ($fasta_infile, $fasta_header_list_file, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'l=s'    => \$fasta_header_list_file,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $fasta_header_list_file
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -l fasta_header_list_file -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 
    
      -l fasta_header_list_file - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_filename = fileparse($fasta_infile);
my $fasta_header_list_filename = fileparse($fasta_header_list_file);

warn "Loading $fasta_header_list_filename ids for parsing $fasta_filename ....\n\n";
my $header_id_list = parse_header_id_list($fasta_header_list_file);

my %fasta_seqs = ();
my $fasta_outfile = join('/', $output_dir, join("_", $fasta_filename, $fasta_header_list_filename . ".fasta"));
open(OUTFILE1, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
my $desc_outfile = join('/', $output_dir, join("_", $fasta_filename, $fasta_header_list_filename, "descriptions.txt"));
open(OUTFILE2, ">$desc_outfile") or die "Couldn't open file $desc_outfile for writting, $!"; 
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $seq_desc = $seq_entry->desc;

      my $fasta_unique_id = join(" ", $seq_id, $seq_desc);
#       warn $fasta_unique_id . "\n";
      if(defined($header_id_list->{$fasta_unique_id})){
	    warn join("\n", ">$fasta_unique_id", $sequence) . "\n";
	    print OUTFILE1 join("\n", ">$fasta_unique_id", $sequence) . "\n";
	    print OUTFILE2 "$fasta_unique_id" . "\n";
      }
     
}
close(OUTFILE1) or die "Couldn't close file $fasta_outfile"; 
close(OUTFILE2) or die "Couldn't close file $desc_outfile";

# Clean out the sequence I/O object.
$seqio = ();

sub parse_header_id_list{
      my $fasta_header_list_file = shift;
      die "Error lost organism ids input file" unless defined $fasta_header_list_file;

      my %organism_list;
      open(INFILE, "<$fasta_header_list_file") or die "Couldn't open file $fasta_header_list_file for reading, $!";
      while(<INFILE>){
	    chomp $_;
	    warn "$_\n";
	    my $organism_id = $_;
	    $organism_list{$organism_id} = $organism_id;

      }
      close(INFILE) or die "Couldn't close file $fasta_header_list_file";
      
      return \%organism_list;
}