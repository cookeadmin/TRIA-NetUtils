#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
# perl get_fasta_seqs.pl -i /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_db_2014-05-30 -o /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_organisms_db

my ($fasta_infile, $organism_ids_infile, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}




my $fasta_filename = fileparse($fasta_infile);


my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
my $fasta_outfile = join('/', $output_dir, $fasta_filename);
open(OUTFILE, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = uc($seq_entry->seq);
      my $seq_desc = $seq_entry->desc;
      my $sequence_length = length($sequence);

      my $fasta_unique_id = join(" ", $seq_id, $seq_desc, "length=$sequence_length");
      warn $fasta_unique_id . "\n";
      my $trimmed_sequence = get_subseq($sequence, 1, 33);
      print OUTFILE join("\n", join("", ">", $fasta_unique_id), $trimmed_sequence) . "\n";
      
}
close(OUTFILE) or die "Couldn't close file $fasta_outfile"; 
# Clean out the sequence I/O object.
$seqio = ();

# my $seq = get_subseq("AGCTTGCGTT", 3, 8);
# warn $seq . "\n";
sub get_subseq{

        my $sequence = shift;
        die "Error lost sequence" unless defined $sequence;

        my $seq_start = shift;
        die "Error lost start of sequence" unless defined $seq_start;

        my $seq_end = shift;
        die "Error lost end of sequence" unless defined $seq_end;

        $seq_start = $seq_start - 1;
        $seq_end = $seq_end;

        my $length = ($seq_end - $seq_start);

        my $trimmed_seq = substr($sequence, $seq_start, $length);

        return uc($trimmed_seq);
}
