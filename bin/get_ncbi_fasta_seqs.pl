#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
# perl get_fasta_seqs.pl -i /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_db_2014-05-30 -l ~/workspace/adriana/organisms-list.txt -o /home/tria-assembly-archive/blast_databases/ncbi_nr_db/ncbi_nr_organisms_db

my ($fasta_infile, $organism_ids_infile, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'l=s'    => \$organism_ids_infile,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $organism_ids_infile
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -l organism_ids_infile -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 
    
      -l organism_ids_infile - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}




my $fasta_filename = fileparse($fasta_infile);
my $organism_ids_filename = fileparse($organism_ids_infile);

warn "Loading $organism_ids_filename organism name-ids for parsing $fasta_filename ....\n\n";
my $organism_list = parse_organism_list($organism_ids_infile);


my $fasta_outfile = join('/', $output_dir, join("_", $fasta_filename, $organism_ids_filename . ".fasta"));
open(OUTFILE1, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
my $desc_outfile = join('/', $output_dir, join("_", $fasta_filename, $organism_ids_filename, "descriptions.txt"));
open(OUTFILE2, ">$desc_outfile") or die "Couldn't open file $desc_outfile for writting, $!";
my %fasta_seqs = ();
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $seq_desc = $seq_entry->desc;


      my @split_seq_desc = split(/\001/, $seq_desc);
      my %organism_ids = ();
      foreach my $seq_dec_entry (@split_seq_desc){
	    if($seq_dec_entry =~ m/\[(.+)\]/){
		  my $organism_id = $1;
		  $organism_ids{$organism_id} = $organism_id;
	    }
      }

      my $fasta_unique_id = join(" ", $seq_id, $seq_desc);
      foreach my $organism_id (sort keys %organism_ids){
	    if(defined($organism_list->{$organism_id})){
		  warn $organism_id . "\n";
# 			warn $seq_id . "\n";
# 			warn $seq_entry->desc . "\n";
		  warn $fasta_unique_id . "\n";
		  warn $sequence . "\n";
		  
		  print OUTFILE1 join("\n", ">$fasta_unique_id", $sequence) . "\n";
		  print OUTFILE2 "$fasta_unique_id" . "\n";
		  last;
	    }
      }
}
close(OUTFILE1) or die "Couldn't close file $fasta_outfile"; 
close(OUTFILE2) or die "Couldn't close file $desc_outfile";

# Clean out the sequence I/O object.
$seqio = ();

sub parse_organism_list{
      my $organism_ids_infile = shift;
      die "Error lost organism ids input file" unless defined $organism_ids_infile;

      my %organism_list;
      open(INFILE, "<$organism_ids_infile") or die "Couldn't open file $organism_ids_infile for reading, $!";
      while(<INFILE>){
	    chomp $_;
# 	    warn "$_\n";
	    my $organism_id = $_;
	    $organism_list{$organism_id} = $organism_id;

      }
      close(INFILE) or die "Couldn't close file $organism_ids_infile";
      
      return \%organism_list;
}