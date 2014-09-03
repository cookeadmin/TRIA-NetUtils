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

my %fasta_seqs;
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = uc($seq_entry->seq);
      my $seq_desc = $seq_entry->desc;
      my $sequence_length = length($sequence);

      my $fasta_unique_id = join(" ", $seq_id, $seq_desc, "length=$sequence_length");
      warn $fasta_unique_id . "\n";
      $sequence =~ s/X/N/g;
      warn $sequence . "\n";
      $fasta_seqs{$fasta_unique_id} = $sequence;
}

# Clean out the sequence I/O object.
$seqio = ();

my $fasta_outfile = join('/', $output_dir, $fasta_filename);
open(OUTFILE, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
foreach my $fasta_unique_id (sort keys %fasta_seqs){
      print OUTFILE ">$fasta_unique_id" . "\n";
      print OUTFILE $fasta_seqs{$fasta_unique_id} . "\n";
}
close(OUTFILE) or die "Couldn't close file $fasta_outfile"; 
