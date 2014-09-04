#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use Bio::SearchIO;

my ($target_infile, $query_infile, $min_percent_id, $num_descriptions, $num_alignments, $blast_num_cpu, $output_dir);
GetOptions(
      'd=s'    => \$target_infile,
      'i=s'    => \$query_infile,
      'p=s'    => \$min_percent_id,
      'v=s'    => \$num_descriptions,
      'b=s'    => \$num_alignments,
      'a=s'    => \$blast_num_cpu,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $target_infile
      and defined $query_infile
      and defined $output_dir
);


$min_percent_id = 80 unless defined $min_percent_id;
$num_descriptions = 5 unless defined $num_descriptions;
$num_alignments = 5 unless defined $num_alignments;
$blast_num_cpu = 2 unless defined $blast_num_cpu;

my ($makeblastdb, $tblastx);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$tblastx			= '/usr/local/bin/tblastx';

sub usage {

die <<"USAGE";

Usage: $0 -d target_infile -i query_infile -p min_percent_id -v num_descriptions -b num_alignments -a blast_num_cpu -o output_dir

Description - 

OPTIONS:

      -d target_infile -

      -i query_infile -

      -p min_percent_id 

      -v num_descriptions -

      -b num_alignments -

      -a blast_num_cpu -

      -o output_dir -


USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_query_name = fileparse($query_infile);
my $fasta_target_name = fileparse($target_infile);
my $tblastx_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . '.tblastx'));

my $tblastx_infile = generate_tblastx($query_infile, $target_infile, $num_descriptions, $num_alignments, $min_percent_id, $blast_num_cpu, $tblastx_filename);

my %query_annotations = ();

open(INFILE, "<$query_infile") or die "Couldn't open file $query_infile for reading, $!";
while(<INFILE>){
      chomp $_;

      if($_ =~ m/^>(.+)/){
	    my $query_header = $1;
# 	    warn $query_header . "\n";
	    my ($query_id, $query_desc) = split(/\s{1}|\001/, $query_header, 2);
# 	    warn join(" ===> ", $query_id, $query_header) . "\n";
	    $query_annotations{$query_id} = $query_header;
      }
}
close(INFILE) or die "Couldn't close file $query_infile";

my %target_annotations;
open(INFILE, "<$target_infile") or die "Couldn't open file $target_infile for reading, $!";
while(<INFILE>){
      chomp $_;

      if($_ =~ m/^>(.+)/){
	    my $target_header = $1;
# 	    warn $target_header . "\n";
	    my ($target_id, $target_desc) = split(/\s{1}|\001/, $target_header, 2);
# 	    warn join(" ===> ", $target_id, $target_header) . "\n";
	    $target_annotations{$target_id} = $target_header;
      }
}
close(INFILE) or die "Couldn't close file $target_infile";

warn "Generating tblastx tsv file....\n";
my $tblastx_outfile = $tblastx_filename . ".tsv";
open(OUTFILE, ">$tblastx_outfile") or die "Couldn't open file $tblastx_outfile for writting, $!";
print OUTFILE join("\t", "query_name", "target_name", "percent_identity", "align_length", "num_mismatch", 
      "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n"; 
my $searchio = Bio::SearchIO->new(-file   => $tblastx_infile,
                                  -format => 'blast') or die "Error: Can't parse blast file $tblastx_infile using SearchIO $!";
while(my $result = $searchio->next_result){
      while(my $hit = $result->next_hit){
	    while(my $hsp = $hit->next_hsp){
		  if ($hsp->percent_identity >= $min_percent_id) {
			my ($query_name,$target_name,$percent_identity,$align_length,$num_mismatch,$num_gaps,$query_strand,
			      $query_start,$query_end,$target_start,$target_end,$target_strand,$e_value,$bit_score);
			$query_name = $query_annotations{$result->query_name};
			$target_name = $target_annotations{$hit->name};
			$percent_identity = sprintf("%.2f", $hsp->percent_identity);
			$align_length= $hsp->length('total');
			$num_mismatch = ($hsp->length('total') - ($hsp->num_identical + $hsp->gaps));
			$num_gaps = $hsp->gaps;

			$query_strand = $hsp->strand('query');
			if($query_strand > 0){
			      $query_start = $hsp->start('query');
			      $query_end = $hsp->end('query');
			}elsif($query_strand < 0){
			      $query_start = $hsp->end('query');
			      $query_end = $hsp->start('query');
			}
			
			$target_strand = $hsp->strand('hit');
			if($target_strand > 0){
			      $target_start = $hsp->start('hit');
			      $target_end = $hsp->end('hit');
			}elsif($target_strand < 0){
			      $target_start = $hsp->end('hit');
			      $target_end = $hsp->start('hit');
			}

			$e_value = $hsp->evalue;
			$e_value = "< 1e-179" if ($e_value =~ m/0\.0/);
			$bit_score = $hsp->bits;
			print OUTFILE join("\t", $query_name, $target_name, $percent_identity, $align_length, $num_mismatch, 
			      $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
		  }
	    }
      }
}
close(OUTFILE) or die "Couldn't close file $tblastx_outfile";

# makeblastdb -in ncbi_nr_db_2014-05-30_organisms.fasta -dbtype 'nucl' -out ncbi_nr_db_2014-05-30_organisms.fasta
sub makeblastdb_nuc{
      my $fastadb = shift;
      die "Error lost fastadb to makeblastdb" unless defined $fastadb;

      # format the database file into .nin .nsq .nhr files.
      my ($fastadbNIN, $fastadbNSQ, $fastadbNHR);
      $fastadbNIN = $fastadb . '.nin';
      $fastadbNSQ = $fastadb . '.nsq';
      $fastadbNHR = $fastadb . '.nhr';
      unless(-s $fastadbNIN and -s $fastadbNSQ and -s $fastadbNHR){
	    warn "Calling makeblastdb for $fastadb....\n";
	    warn "$makeblastdb -in $fastadb -dbtype nucl\n\n";
	    system($makeblastdb, 
		  '-in', $fastadb, 
		  '-dbtype', 'nucl'
	    ) == 0 or die "Error calling $makeblastdb -in $fastadb -dbtype nucl: $?";
      }

}

sub generate_tblastx{
      my $fasta_query = shift;
      die "Error lost fasta query file" unless defined $fasta_query;
      my $fasta_target = shift;
      die "Error lost fasta database target file" unless defined $fasta_target;
      my $num_descriptions = shift;
      die "Error lost number of descriptions" unless defined $num_descriptions;
      my $num_alignments = shift;
      die "Error lost number of alignments" unless defined $num_alignments;
      my $min_percent_id = shift;
      die "Error lost minimum percent identity" unless defined $min_percent_id;
      my $blast_num_cpu = shift;
      die "Error lost number of cpus to allocate" unless defined $blast_num_cpu;
      my $tblastx_filename = shift;
      die "Error lost tblastx output filename" unless defined $tblastx_filename;

      makeblastdb_nuc($fasta_target);

      my $tblastx_outfile = $tblastx_filename . ".aln";
      unless(-s $tblastx_outfile){
	    warn "Generating tblastx file....\n";
	    my $tblastxCmd  = "$tblastx -query $fasta_query -db $fasta_target -dust yes -num_descriptions $num_descriptions -num_alignments $num_alignments -evalue 1e-6 -out $tblastx_outfile -num_threads $blast_num_cpu";
	    warn $tblastxCmd . "\n\n";

	    my $status = system($tblastx, 
		  '-query', $fasta_query,
		  '-db', $fasta_target,
		  '-dust', 'yes',
		  '-num_descriptions', $num_descriptions,
		  '-num_alignments', $num_alignments,
		  '-evalue', 1e-6,
		  '-out', $tblastx_outfile,
		  '-num_threads', $blast_num_cpu
	    ) == 0 or die "Error calling $tblastx: $?";

      }
      return $tblastx_outfile;
}
