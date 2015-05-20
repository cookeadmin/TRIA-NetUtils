#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use List::MoreUtils qw( minmax );
use POSIX qw/strftime/;
use File::Basename;
use Bio::SeqIO;

my ($splign_infile, $refgen_infile, $output_dir);
GetOptions(
      'i=s'    => \$splign_infile,
      'g=s'    => \$refgen_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $splign_infile
      and defined $refgen_infile
      and defined $output_dir
);

sub usage {

die <<"USAGE";

Usage: $0 -i splign_infile -g refgen_infile -o output_dir

Description - 

OPTIONS:

      -i splign_infile -
      
      -g refgen_infile -
      
      -o output_dir -

USAGE
}

my $date = strftime('%Y-%m-%d', localtime);

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $genome_filename = fileparse($refgen_infile, qr/\.fasta/);
my %fasta_seqs = ();
my $seqio = Bio::SeqIO->new(-file => $refgen_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $seq_desc = $seq_entry->desc;

      my $fasta_unique_id = join(" ", $seq_id, $seq_desc);

      my @split_fasta_unique_id = split(/_/, $fasta_unique_id);
      my $target_id = $split_fasta_unique_id[0];
      
      $fasta_seqs{$target_id} = join("\n", ">$fasta_unique_id", $sequence);
}

my %splign_entries = ();
open(INFILE,"<$splign_infile") or die "Error opening $splign_infile: $!";
while(<INFILE>){
	chomp $_;
# 	warn $_ . "\n";
	my @split_splign_entry = split(/\t/, $_);
	
	my ($query_strand_id, $query_id, $target_id, $percent_id, $align_length, $query_start, $query_end, $target_start, $target_end, $align_type, $align_transcript) = @split_splign_entry;
	my $query_strand = $query_strand_id;
	$query_strand =~ s/\d+//g;
	
	my $splign_align_id = join("\001", $query_id, $target_id, $query_strand_id);
	my $parsed_splign_entry = join("\t", $query_strand, $query_id, $target_id, $percent_id, $align_length, $query_start, $query_end, $target_start, $target_end, $align_type, $align_transcript);
	push(@{$splign_entries{$splign_align_id}}, $parsed_splign_entry);
	
}
close INFILE or die "Error closing $splign_infile: $!";


my %filtered_splign = ();
foreach my $splign_align_id (sort {$a cmp $b} keys %splign_entries){

# 	warn $splign_align_id . "\n";
	my ($aligned_counts, $unaligned_counts) = 0;
	my (@query_coords, @target_coords) = ();
	foreach my $splign_entry (@{$splign_entries{$splign_align_id}}){
		my @split_splign_entry = split(/\t/, $splign_entry);
		my ($query_strand, $query_id, $target_id, $percent_id, $align_length, $query_start, $query_end, $target_start, $target_end, $align_type, $align_transcript) = @split_splign_entry;
		
		if($align_type =~ m/<exon>/){
			my $alignment_length = (abs($query_end - $query_start) + 1);
			$aligned_counts += $alignment_length;
		}
		
		if($align_type !~ m/<exon>/){
			my $alignment_length = (abs($query_end - $query_start) + 1);
			$unaligned_counts += $alignment_length;
		}
		
		push(@query_coords, $query_start, $query_end) if($align_type =~ m/<exon>/);
		push(@target_coords, $target_start, $target_end) if($align_type =~ m/<exon>/);
		print $splign_entry . "\n";
	}
	
	$unaligned_counts = 0 unless(defined($unaligned_counts));
	$aligned_counts = 0 unless(defined($aligned_counts));
	
	
	my ($query_start, $query_end) = minmax(@query_coords);
	my ($target_start, $target_end) = minmax(@target_coords);
	
	$filtered_splign{$splign_align_id} = join("\t", $aligned_counts, $unaligned_counts, $query_start, $query_end, $target_start, $target_end);
}

my $splign_filename = fileparse($splign_infile, qr/\.txt/);
my $splign_outfile = join('/', $output_dir, $splign_filename . ".fasta");
open(OUTFILE1, ">$splign_outfile") or die "Couldn't open file $splign_outfile for writting, $!";
my $splign_upstream_outfile = join('/', $output_dir, join("-", $splign_filename, "upstream.fasta"));
open(OUTFILE2, ">$splign_upstream_outfile") or die "Couldn't open file $splign_upstream_outfile for writting, $!";
my $splign_downstream_outfile = join('/', $output_dir, join("-", $splign_filename, "downstream.fasta"));
open(OUTFILE3, ">$splign_downstream_outfile") or die "Couldn't open file $splign_downstream_outfile for writting, $!";
foreach my $splign_align_id (sort {$a cmp $b} keys %filtered_splign){

# 	warn join("\t", $splign_align_id, $filtered_splign{$splign_align_id}) . "\n";
	

# 		my $splign_align_id = join("\001", $query_id, $target_id, $query_strand_id);
	my @split_splign_align_id = split(/\001/, $splign_align_id);
	my ($query_id, $target_id, $query_strand_id) = @split_splign_align_id;
	
	my @split_splign_entry = split(/\t/, $filtered_splign{$splign_align_id});
	my ($aligned_counts, $unaligned_counts, $query_start, $query_end, $target_start, $target_end) = @split_splign_entry;
	
	my @split_target_unique_id = split(/_/, $target_id);
	my $target_name = $split_target_unique_id[0];
	
	my $genomic_sequence_entry = $fasta_seqs{$target_name};
	my ($genomic_seq_id, $genomic_sequence) = split(/\n/, $genomic_sequence_entry, 2);

# 			next if ($genomic_sequence =~ m/N/);
	
# 			get_subseq($genomic_sequence, 1, ($target_start - 1));
	my $genomic_cDNA_sequence = get_subseq($genomic_sequence, $target_start, $target_end);
	
	$genomic_cDNA_sequence =~ s/(.{80})/$1\n/gs;

	my $mRNA_length = (($target_end - $target_start) + 1);
	my $cDNA_length = ($aligned_counts + $unaligned_counts);
	
	my $q_cov = sprintf("%.2f", ($aligned_counts/($aligned_counts + $unaligned_counts)) * 100);
	my $query_coverage =  $q_cov . "%";
	
	my $strand = $query_strand_id;
	$strand =~ s/\d+//g;
	
	print OUTFILE1 join(" ", join(":", ">$query_id", "Splign", $genome_filename), "type=genomic", "query_coverage=$query_coverage", join(":", "gene=$query_id", $query_start, $query_end, $strand, "length=$cDNA_length"), join(":", "source=contig$target_name", $target_start, $target_end, "+", "length=$mRNA_length"), "date=$date") . "\n";
	print OUTFILE1 $genomic_cDNA_sequence . "\n\n";
	
	my $genomic_upstream_target_start = 1;
	my $genomic_upstream_target_end = ($target_start - 1);
	
	my $genomic_upstream_sequence = get_subseq($genomic_sequence, $genomic_upstream_target_start, $genomic_upstream_target_end);
	
	my $genomic_upstream_sequence_length = length($genomic_upstream_sequence);
	my ($genomic_upstream_seq, $genomic_upstream_seq_start, $genomic_upstream_seq_end, $genomic_upstream_seq_length);
	if($genomic_upstream_sequence_length > 4000){
	
		($genomic_upstream_seq_start, $genomic_upstream_seq_end) = ((($genomic_upstream_target_end - 4000) + 1), $genomic_upstream_target_end);
		$genomic_upstream_seq = get_subseq($genomic_upstream_sequence, (($genomic_upstream_target_end - 4000) + 1), $genomic_upstream_target_end);
		$genomic_upstream_seq_length = length($genomic_upstream_seq);
	}elsif($genomic_upstream_sequence_length <= 4000){
		($genomic_upstream_seq_start, $genomic_upstream_seq_end) = ($genomic_upstream_target_start, $genomic_upstream_target_end);
		$genomic_upstream_seq = $genomic_upstream_sequence;
		$genomic_upstream_seq_length = length($genomic_upstream_seq);
	}
	
	if($genomic_upstream_seq_length ne 0){
		$genomic_upstream_seq =~ s/(.{80})/$1\n/gs;
		
		if($strand eq "+"){
			print OUTFILE2 join(" ", join(":", ">$query_id", $strand, "Splign", $genome_filename), "type=genomic upstream sequence", join(":", "source=contig$target_name", $genomic_upstream_seq_start, $genomic_upstream_seq_end, "+", "length=$genomic_upstream_seq_length"), "date=$date") . "\n";
			print OUTFILE2 $genomic_upstream_seq . "\n\n";
		}elsif($strand eq "-"){
			
			print OUTFILE3 join(" ", join(":", ">$query_id", $strand, "Splign", $genome_filename), "type=genomic downstream sequence", join(":", "source=contig$target_name", $genomic_upstream_seq_start, $genomic_upstream_seq_end, "+", "length=$genomic_upstream_seq_length"), "date=$date") . "\n";
			print OUTFILE3 $genomic_upstream_seq . "\n\n";
		}
	}
	
	my $genomic_downstream_target_start = ($target_end + 1);
	my $genomic_downstream_target_end = length($genomic_sequence);
	my $genomic_downstream_sequence = get_subseq($genomic_sequence, $genomic_downstream_target_start, $genomic_downstream_target_end);
	
	my $genomic_downstream_sequence_length = length($genomic_downstream_sequence);
	
	my ($genomic_downstream_seq, $genomic_downstream_seq_start, $genomic_downstream_seq_end, $genomic_downstream_seq_length);
	if($genomic_downstream_sequence_length > 1000){
		($genomic_downstream_seq_start, $genomic_downstream_seq_end) = ($genomic_downstream_target_start, (($genomic_downstream_target_start + 1000) - 1));
		$genomic_downstream_seq = get_subseq($genomic_sequence, $genomic_downstream_target_start, (($genomic_downstream_target_start + 1000) - 1));
		$genomic_downstream_seq_length = length($genomic_downstream_seq);
	}elsif($genomic_downstream_sequence_length <= 1000){
		($genomic_downstream_seq_start, $genomic_downstream_seq_end) = ($genomic_downstream_target_start, $genomic_downstream_target_end);
		$genomic_downstream_seq = get_subseq($genomic_sequence, $genomic_downstream_target_start, $genomic_downstream_target_end);
		$genomic_downstream_seq_length = length($genomic_downstream_seq);
	}
	
	if($genomic_downstream_seq_length ne 0){
		$genomic_downstream_seq =~ s/(.{80})/$1\n/gs;
		
		if($strand eq "+"){
			print OUTFILE3 join(" ", join(":", ">$query_id", $strand, "Splign", $genome_filename), "type=genomic downstream sequence", join(":", "source=contig$target_name", $genomic_downstream_seq_start, $genomic_downstream_seq_end, "+", "length=$genomic_downstream_seq_length"), "date=$date") . "\n";
			print OUTFILE3 $genomic_downstream_seq . "\n\n";
		}elsif($strand eq "-"){
			print OUTFILE2 join(" ", join(":", ">$query_id", $strand, "Splign", $genome_filename), "type=genomic upstream sequence", join(":", "source=contig$target_name", $genomic_downstream_seq_start, $genomic_downstream_seq_end, "+", "length=$genomic_downstream_seq_length"), "date=$date") . "\n";
			print OUTFILE2 $genomic_downstream_seq . "\n\n";
		}
	}
}
close(OUTFILE1) or die "Couldn't close file $splign_outfile";
close(OUTFILE2) or die "Couldn't close file $splign_upstream_outfile";
close(OUTFILE3) or die "Couldn't close file $splign_downstream_outfile";

# $trimmed_seq = get_subseq($sequence, $seq_start, $seq_end) - Get the subsequence based on the input sequence, sequence start, and sequence end.
#
# Input paramater(s):
#
# $sequence - The input sequence to obtain a subsequence.
#
# $seq_start - The start position of the sequence.
#
# $seq_end - The end position of the sequence.
#
# Output paramater(s):
#
# $trimmed_seq - The trimmed subsequence of the input sequence.
sub get_subseq{

	# The input sequence to obtain a subsequence.
        my $sequence = shift;
        die "Error lost input sequence to obtain a subsequence" unless defined $sequence;

        # The start position of the sequence.
        my $seq_start = shift;
        die "Error lost start position of the sequence" unless defined $seq_start;

        # The end position of the sequence.
        my $seq_end = shift;
        die "Error lost end position of the sequence" unless defined $seq_end;

        $seq_start = $seq_start - 1;
        $seq_end = $seq_end;

        my $length = ($seq_end - $seq_start);

        my $trimmed_seq = substr($sequence, $seq_start, $length);

        return uc($trimmed_seq);
}

