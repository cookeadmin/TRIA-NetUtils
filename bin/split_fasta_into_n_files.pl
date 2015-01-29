#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
# perl get_fasta_seqs.pl -i /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_db_2014-05-30 -n 10 -o /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_organisms_db

my ($fasta_infile, $num_files, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'n=s'    => \$num_files,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $fasta_infile
      and defined $output_dir
);

$num_files = 10 unless defined $num_files;

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -n num_files -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 
    
      -n num_files - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_file_parts = split_fasta_into_n_files($fasta_infile, $num_files, $output_dir);



sub split_fasta_into_n_files{

      my $fasta_infile = shift;
      die "Error lost fasta input file" unless defined $fasta_infile;

      my $num_files = shift;
      die "Error lost organism ids input file" unless defined $num_files;

      my $output_dir = shift;
      die "Error lost output directory" unless defined $output_dir;
      
      my $fasta_filename = fileparse($fasta_infile, qr/\.\w+/);

      warn "Splitting $fasta_filename into $num_files files....\n\n";

      my @fasta_seqs;
      my $num_sequences = 0;
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
	    push(@fasta_seqs, $fasta_seq);
	    $num_sequences++;
      }
     
      # Clean out the sequence I/O object.
      $seqio = ();

      my $num_seq_per_file = int($num_sequences / $num_files);
      my $remainder_seq_per_file = int($num_sequences % $num_files);
      my %fasta_file_parts;
      my $seq_counter = 0;
      for(my $i = 1; $i <= $num_files; $i++){
	    
	    my $num_seqs;
	    if(($i >= 1) and ($i < $num_files)){

		  $num_seqs = ($num_seq_per_file * $i);

	    }else{
		  
		  $num_seqs = (($num_seq_per_file * $i) + $remainder_seq_per_file);
	    }

	    warn join(" ", "File $i:", join("-", $seq_counter, ($num_seqs - 1))) . "\n";
	    my $fasta_file = join("-", $fasta_filename, "part$i");
	    my $fasta_outfile = join('/', $output_dir, $fasta_file . ".fasta");

	    $fasta_file_parts{$fasta_file} = $fasta_outfile;

	    open(OUTFILE, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
	    for(my $j = $seq_counter; $j < $num_seqs; $j++){
		  print OUTFILE $fasta_seqs[$j]. "\n";
	    }
	    close(OUTFILE) or die "Couldn't close file $fasta_outfile";
	    $seq_counter = $num_seqs;
      }
      return \%fasta_file_parts;
}