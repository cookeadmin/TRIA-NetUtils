#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
# perl filter_genome_seqs.pl -i /home/tria-assembly-archive/white_spruce/WS77111-scaffolds.fa -l 2000 -o ~/workspace/WS77111-scaffolds 
my ($fasta_infile, $min_length, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'l=s'    => \$min_length,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $output_dir
);

$min_length = 2000 unless(defined($min_length));

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -l min_length -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 
    
      -l min_length - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_filename = fileparse($fasta_infile);

my %fasta_seqs = ();
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
my $fasta_outfile = join('/', $output_dir, join("_", $fasta_filename, join("-", "min", "length", join("", $min_length, "bp")) . ".fasta"));
open(OUTFILE1, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
my $reads_outfile = join('/', $output_dir, join("_", $fasta_filename, "all-reads.info"));
open(OUTFILE2, ">$reads_outfile") or die "Couldn't open file $reads_outfile for writting, $!";
my $min_reads_outfile = join('/', $output_dir, join("_", $fasta_filename, join("-", "min", "length", join("", $min_length, "bp"), "reads") . ".info"));
open(OUTFILE3, ">$min_reads_outfile") or die "Couldn't open file $min_reads_outfile for writting, $!";
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $seq_desc = $seq_entry->desc;

      my $fasta_unique_id = join(" ", $seq_id, $seq_desc);

      my $sequence_length = length($sequence);

      my $sequence_id_length = join("\t", $fasta_unique_id, $sequence_length);
      print OUTFILE2 $sequence_id_length . "\n";
      if($sequence_length >= $min_length){
# 	    warn $fasta_unique_id . "\n";
# 	    warn $sequence . "\n";
	    print OUTFILE1 join("\n", ">$fasta_unique_id", $sequence) . "\n";
	    print OUTFILE3 $sequence_id_length . "\n";
      }

}
close(OUTFILE1) or die "Couldn't close file $fasta_outfile"; 
close(OUTFILE2) or die "Couldn't close file $reads_outfile";
close(OUTFILE3) or die "Couldn't close file $min_reads_outfile";

# Clean out the sequence I/O object.
$seqio = ();
