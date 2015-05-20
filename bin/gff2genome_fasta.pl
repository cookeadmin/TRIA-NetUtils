#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use POSIX qw/strftime/;
use File::Basename;
use Bio::SeqIO;

my ($refgen_infile, $gff_infile, $output_dir);
GetOptions(
      'i=s'    => \$refgen_infile,
      'g=s'    => \$gff_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $refgen_infile
      and $gff_infile
      and defined $output_dir
);

sub usage {

die <<"USAGE";

Usage: $0 -i refgen_infile -g gff_infile -o output_dir

Description - 

OPTIONS:

      -i refgen_infile -
      -g gff_infile - 
      -o output_dir -

USAGE
}

my $gff_filename = fileparse($gff_infile, qr/.gff/);

my ($genome_filename, $cDNA_filename) = split(/_/, $gff_filename, 2);

my $date = strftime('%Y-%m-%d', localtime);


# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

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

my %gff_hits = ();
my $i = 0;
open(INFILE,"<$gff_infile") or die "Error opening $gff_infile: $!";
while(<INFILE>){
	chomp $_;
	if($i ne 0){
# 		warn $_ . "\n";
		my @split_blastx_entry = split(/\t/, $_);
		
		# 142155638       BLAT    cDNA_match      182     2326    98.48   -       .       ID=GQ03808_D10.1_mid2;Target=GQ03808_D10.1 1121 1
		my ($target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_blastx_entry;
		if($match_type eq "cDNA_match"){
			my $gff_entry = join("\t", $target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
			
			my $query_id = "";
			if($attributes =~ m/ID=(.+);Target=.+\s\d+\s\d+/){ # ID=WS00720_I24.1_mid7;Target=WS00720_I24.1 908 30
				$query_id = $1;
			}else{
				die "Error: Attributes not in correct format: $attributes";
			}
			
			my $unique_id = join("_", $target_id, $query_id);
			$gff_hits{$unique_id}{"cDNA_match"} = $gff_entry;
		}
		if($match_type =~ m/exon/){
			my $gff_entry = join("\t", $target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
			
			my $query_id = "";
			if($attributes =~ m/Parent=(.+);Target=.+\s\d+\s\d+/){ # Parent=WS0337_D11.1_mid1;Target=WS0337_D11.1 20 6
				$query_id = $1;
			}else{
				die "Error: Attributes not in correct format: $attributes";
			}
			
			my $unique_id = join("_", $target_id, $query_id);
			push(@{$gff_hits{$unique_id}{"exon"}}, $gff_entry);
		}
	}
	$i++;
}
close INFILE or die "Error closing $gff_infile: $!";


my %filtered_gff = ();
foreach my $unique_id (sort {$a cmp $b} keys %gff_hits){
	my @split_blastx_entry = split(/\t/, $gff_hits{$unique_id}{"cDNA_match"});
	my ($target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_blastx_entry;

	my $query_id;
	if($attributes =~ m/ID=.+;Target=(.+)\s\d+\s\d+/){ # ID=WS00720_I24.1_mid7;Target=WS00720_I24.1 908 30
		$query_id = $1;
	}else{
		die "Error: Attributes not in correct format: $attributes";
	}
	my @gff_entries = ();
	push(@gff_entries, $gff_hits{$unique_id}{"cDNA_match"});

	foreach my $exon_entry (@{$gff_hits{$unique_id}{"exon"}}){
	
		my @split_blastx_entry = split(/\t/, $exon_entry);
		my ($target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_blastx_entry;
		
		my $new_exon_entry = join("\t", $target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
		push(@gff_entries, $new_exon_entry);

	}
	
	my $new_gff_entry = join("\n", @gff_entries);
	push(@{$filtered_gff{$query_id}{$percent_id}}, $new_gff_entry);
	
}


# my $blat2gff_outfile = join('/', $output_dir, join("-", $gff3_filename, "contig-sorted.gff"));
# open(OUTFILE, ">$blat2gff_outfile") or die "Couldn't open file $blat2gff_outfile for writting, $!";
my $counter = 0;
foreach my $query_id (sort {$a cmp $b} keys %filtered_gff){
	foreach my $percent_id (sort {$b <=> $a} keys $filtered_gff{$query_id}){
		foreach my $gff_entry (@{$filtered_gff{$query_id}{$percent_id}}){
			my @split_gff_entry = split(/\n/, $gff_entry);
			my $cDNA_match = $split_gff_entry[0];
			my @split_cDNA_match = split(/\t/, $cDNA_match);
			my ($target_name, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_cDNA_match;
			
#  			die join("\t", $target_name, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
			my ($query_name, $query_start, $query_end);
			if($attributes =~ m/ID=.+;Target=(.+)\s(\d+)\s(\d+)/){ # ID=WS00720_I24.1_mid7;Target=WS00720_I24.1 908 30
				($query_name, $query_start, $query_end) = ($1, $2, $3) if($strand eq "+");
				($query_name, $query_start, $query_end) = ($1, $2, $3) if($strand eq ".");
				($query_name, $query_start, $query_end) = ($1, $3, $2) if($strand eq "-");
				
			}else{
				die "Error: Attributes not in correct format: $attributes";
			}
# 	

			my $genomic_sequence_entry = $fasta_seqs{$target_name};
			my ($genomic_seq_id, $genomic_sequence) = split(/\n/, $genomic_sequence_entry, 2);
			
# 			next if ($genomic_sequence =~ m/N/);
			
# 			get_subseq($genomic_sequence, 1, ($target_start - 1));
			my $genomic_cDNA_sequence = get_subseq($genomic_sequence, $target_start, $target_end);
			$genomic_cDNA_sequence =~ s/(.{80})/$1\n/gs;
# 			get_subseq($genomic_sequence, ($target_end + 1), length($genomic_sequence));

			my $mRNA_length = ($target_end - $target_start);
			my $cDNA_length = ($query_end - $query_start);
			print join(" ", join(":", ">$query_name", $program_source, $genome_filename . ".1"), "type=genomic", "blat_score=$percent_id", join(":", "gene=$query_name", $query_start, $query_end, $strand, "length=$cDNA_length"), join(":", "source=contig$target_name", $target_start, $target_end, "+", "length=$mRNA_length"), "date=$date") . "\n";
			print $genomic_cDNA_sequence . "\n\n";
			# 			GQ03323_E13:genomic.PG29.1 type=genomic gene=GQ03323_E13 source=contig27824:+ date=2013-11-30
			
# 			print OUTFILE $gff_entry . "\n";

# 			last; # We want the best alignments for the fasta files.
		}
# 		last; # We want the best alignments for the fasta files.
	}
	$counter++;
}
# close(OUTFILE) or die "Couldn't close file $blat2gff_outfile";

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

