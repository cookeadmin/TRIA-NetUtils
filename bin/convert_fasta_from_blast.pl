#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;

# perl convert_fasta_from_blast.pl -i ~/workspace/colleen/jack_pine_loblloly_microarray_contigs/jack_pine_illumina_trinity_23Aug2011.fasta_jack-pine-illumina-trinity-23Aug2011-megablast-loblloly-microarray-contigs-list.fasta -m ~/workspace/microarray/jack_pine_microarray_megablast/loblolly_microarray_all_known_sequences.spot-name-indexed.fasta_jack_pine_illumina_trinity_23Aug2011.fasta.megablast.txt -x ~/workspace/colleen/jack_pine_loblloly_microarray_contigs/jack_pine_illumina_trinity_23Aug2011.fasta_jack-pine-illumina-trinity-23Aug2011-megablast-loblloly-microarray-contigs-list.fasta_TAIR10_pep_20101214_updated.blastx.txt -o ~/workspace/colleen
# perl convert_fasta_from_blast.pl -i ~/workspace/loblolly_pine_microarray_analysis/jack_pine_blastx/jack_pine_illumina_trinity_23Aug2011.fasta_jack-pine-illumina-trinity-23Aug2011-megablast-loblloly-microarray-contigs-list.fasta -m ~/workspace/loblolly_pine_microarray_analysis/jack_pine_megablast/loblolly_microarray_all_known_sequences.spot-name-indexed.fasta_jack_pine_illumina_trinity_23Aug2011.fasta.megablast.txt -o ~/workspace/loblolly_pine_microarray_analysis/jack_pine_annotated_sequences

my ($fasta_infile, $megablast_infile, $blastx_infile, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'm=s'    => \$megablast_infile,
#       'x=s'    => \$blastx_infile,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $fasta_infile
      and defined $megablast_infile
#       and defined $blastx_infile
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -i fasta_infile -m megablast_infile -x blastx_infile -o output_dir
    
Description - 
    
OPTIONS:
      -i fasta_infile - 

      -m megablast_infile - 

#       -x blastx_infile - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}


my $fasta_filename = fileparse($fasta_infile);
my $megablast_filename = fileparse($megablast_infile);
# my $blastx_filename = fileparse($blastx_infile);

my %megablast_entries = ();
open(INFILE, "<$megablast_infile") or die "Couldn't open file $megablast_infile for reading, $!";
my $i = 0;
while(<INFILE>){
      chomp $_;
#       warn "$_\n";
      if($i ne 0){
	    my @split_megablast_entry = split(/\t/, $_);
	    my ($query_name, $target_name, $percent_identity, $align_length, $num_mismatch, 
		  $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @split_megablast_entry;
# 	    warn join("\t",$query_name, $target_name, $percent_identity, $align_length, $num_mismatch, 
# 		  $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";	  
	    my $strand = "";
	    if($target_start < $target_end){
		  $strand = "+";
	    }else{
		  $strand = "-";
	    }
	    

	    my @query_header = split(/\s/, $query_name);
	    my ($microarray_id, $microarray_seq_length) = ($query_header[2],$query_header[3]);
	    $microarray_seq_length =~ s/length=/len:/g;

	    $query_name = join("_", $microarray_id, $microarray_seq_length);
	    my $megablast_annotation = join("\001", $target_name, join("; ", $query_name, $percent_identity, $align_length, $strand, $e_value, $bit_score));

# 	    warn $megablast_annotation . "\n";
	    push(@{$megablast_entries{$target_name}}, $megablast_annotation);
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $megablast_infile";

# my %blastx_entries = ();
# open(INFILE, "<$blastx_infile") or die "Couldn't open file $blastx_infile for reading, $!";
# $i = 0;
# while(<INFILE>){
#       chomp $_;
# #       warn "$_\n";
#       if($i ne 0){
# 	    my @split_blastx_entry = split(/\t/, $_);
# 	    my ($query_name, $target_name, $percent_identity, $align_length, $num_mismatch, 
# 		  $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @split_blastx_entry;
# # 	    warn join("\t",$query_name, $target_name, $percent_identity, $align_length, $num_mismatch, 
# # 		  $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";	  
# 	    my $strand = "";
# 	    if($query_start < $query_end){
# 		  $strand = "+";
# 	    }else{
# 		  $strand = "-";
# 	    }
# 	    
# 
# 	    my @target_header = split(/\|/, $target_name);
# 	    
# 	    for(my $i = 0; $i < scalar(@target_header); $i++){
# 		  $target_header[$i] =~ s/^\s//;
# 		  $target_header[$i] =~ s/\s$//;
# 	    }
# 
# 	    my ($annotation_id, $annotation_desc, $annotation_seq_length) = ($target_header[0],$target_header[2],$target_header[3]);
# 	    
# 	    if($annotation_seq_length =~ m/(LENGTH=\d+)/){
# 		  $annotation_seq_length = $1;
# 	    }
# 	    $annotation_seq_length =~ s/LENGTH=/len:/g;
# 
#  	    $target_name = join("_", $annotation_id, $annotation_seq_length,$annotation_desc);
# 	    my $blastx_annotation = join("; ", $target_name, $percent_identity, $align_length, $strand, $e_value, $bit_score);
# 
# # 	    warn $blastx_annotation . "\n";
# 	    push(@{$blastx_entries{$query_name}}, $blastx_annotation);
#       }
#       $i++;
# }
# close(INFILE) or die "Couldn't close file $blastx_infile";

my %fasta_seqs;
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;

      foreach my $megablast_entry (@{$megablast_entries{$seq_id}}){
# 	    if(defined(@{$blastx_entries{$seq_id}})){
# 		  foreach my $blastx_entry (@{$blastx_entries{$seq_id}}){
# 			
# 			print join("\n", join("\001", $megablast_entry, $blastx_entry), $sequence) . "\n";
# 		  }
# 	    }else{
		  print join("\n", $megablast_entry, $sequence) . "\n";
# 	    }
      }
}

# Clean out the sequence I/O object.
$seqio = ();
