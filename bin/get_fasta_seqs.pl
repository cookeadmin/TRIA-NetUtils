#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

# perl get_fasta_seqs.pl -i ~/workspace/beetle_500k_snps/Genomes/Genome_Builds/MPB_SJ072_Genome_2011/SJ072-20110815.fa -l ~/workspace/beetle_500k_snps/Genomes/Genome_Builds/MPB_scaffolds_remove_from2011Genome/MPB_2011Genome_BacteriaContaminatedScaffolds.txt -e yes -o ~/workspace/beetle_500k_snps/Genomes/Genome_Builds/MPB_SJ072_Genome_2011/SJ072-20110815-filtered-seqs 
my ($fasta_infile, $seq_id_list_file, $exclude, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'l=s'    => \$seq_id_list_file,
      'e=s'    => \$exclude,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $seq_id_list_file
      and defined $output_dir
);
$exclude = 'no' unless defined $exclude;

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -l seq_id_list_file -e exclude -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 
    
      -l seq_id_list_file -
      
      -e exclude -
      
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_filename = fileparse($fasta_infile, qr/\.\w+/);
my $seq_id_list_filename = fileparse($seq_id_list_file, qr/\.\w+/);

warn "Loading $seq_id_list_filename ids for parsing $fasta_filename ....\n\n";
my $seq_id_list = parse_seq_id_list($seq_id_list_file);

my %fasta_seqs = ();
my $fasta_outfile = join('/', $output_dir, join("_", $fasta_filename, $seq_id_list_filename . ".fasta"));
open(OUTFILE1, ">$fasta_outfile") or die "Couldn't open file $fasta_outfile for writting, $!";
my $desc_outfile = join('/', $output_dir, join("_", $fasta_filename, $seq_id_list_filename, "descriptions.txt"));
open(OUTFILE2, ">$desc_outfile") or die "Couldn't open file $desc_outfile for writting, $!"; 
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

	my $seq_id = $seq_entry->id;
	my $sequence = $seq_entry->seq;
	my $seq_desc = $seq_entry->desc;
 
 	my $sequence_id = "";
 	if($seq_desc ne ""){
 		$sequence_id = join(" ", $seq_id, $seq_desc);
 	}else{
 		$sequence_id = $seq_id;
 	}
# 	warn $sequence_id . "\n";

	if($exclude eq 'no'){
		my $is_seq_id_found = "false";
		foreach my $seq_id_entry (sort keys %{$seq_id_list}){
		
			if($sequence_id =~ m/$seq_id_entry/){
				$is_seq_id_found = "true";
				last;
			}
		}
		
		if($is_seq_id_found eq "true"){
			warn join("\n", ">$sequence_id", $sequence) . "\n";
			print OUTFILE1 join("\n", ">$sequence_id", $sequence) . "\n";
			print OUTFILE2 "$sequence_id" . "\n";
		}
	}
	
	if($exclude eq 'yes'){
	
		my $is_seq_id_found = "false";
		foreach my $seq_id_entry (sort keys %{$seq_id_list}){
		
			if($sequence_id =~ m/$seq_id_entry/){
				$is_seq_id_found = "true";
				last;
			}
		}
		
		if($is_seq_id_found eq "false"){
			my $seq_length = join("=", "length", length($sequence));
			my $seq_header = join(" ", $sequence_id, $seq_length);
			warn join("\n", ">$seq_header", $sequence) . "\n";
			print OUTFILE1 join("\n", ">$seq_header", $sequence) . "\n";
			print OUTFILE2 "$seq_header" . "\n";
		}
	}
     
}
close(OUTFILE1) or die "Couldn't close file $fasta_outfile"; 
close(OUTFILE2) or die "Couldn't close file $desc_outfile";

# Clean out the sequence I/O object.
$seqio = ();

sub parse_seq_id_list{
	my $seq_id_list_file = shift;
	die "Error lost sequence ids input file" unless defined $seq_id_list_file;

	my %seq_id_list = ();
# 	my $i = 0;
	open(INFILE, "<$seq_id_list_file") or die "Couldn't open file $seq_id_list_file for reading, $!";
	while(<INFILE>){
		chomp $_;
 		warn "$_\n";
		my $seq_id = $_;
		$seq_id_list{$seq_id} = $seq_id;
# 		$i++;
	}
	close(INFILE) or die "Couldn't close file $seq_id_list_file";
	return \%seq_id_list;
}