#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;

my ($target_infile, $query_infile, $min_percent_id, $num_descriptions, $num_alignments, $blast_num_cpu, $output_fmt, $output_dir);
GetOptions(
      'd=s'    => \$target_infile,
      'i=s'    => \$query_infile,
      'p=s'    => \$min_percent_id,
      'v=s'    => \$num_descriptions,
      'b=s'    => \$num_alignments,
      'a=s'    => \$blast_num_cpu,
      'f=s'    => \$output_fmt,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $target_infile
      and defined $query_infile
      and defined $output_dir
);


$min_percent_id = 0 unless defined $min_percent_id;
$num_descriptions = 25 unless defined $num_descriptions;
$num_alignments = 25 unless defined $num_alignments;
$blast_num_cpu = 2 unless defined $blast_num_cpu;
$output_fmt = 'all' unless defined $output_fmt;

my ($makeblastdb, $tblastn);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$tblastn			= '/usr/local/bin/tblastn';

sub usage {

die <<"USAGE";

Usage: $0 -d target_infile -i query_infile -p min_percent_id -v num_descriptions -b num_alignments -a blast_num_cpu -f output_fmt -o output_dir

Description - 

OPTIONS:

      -d target_infile -

      -i query_infile -

      -p min_percent_id 

      -v num_descriptions -

      -b num_alignments -

      -a blast_num_cpu -
      
      -f output_fmt - 
      
      -o output_dir -


USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $fasta_query_name = fileparse($query_infile);
my $fasta_target_name = fileparse($target_infile);
my $tblastn_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . '.tblastn'));

my $tblastn_infile = generate_tblastn($query_infile, $target_infile, $num_descriptions, $num_alignments, $min_percent_id, $blast_num_cpu, $output_fmt, $tblastn_filename);

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

sub generate_tblastn{
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
	my $output_fmt = shift;
	die "Error lost blastn output format" unless defined $output_fmt;
	my $tblastn_filename = shift;
	die "Error lost tblastn output filename" unless defined $tblastn_filename;

	makeblastdb_nuc($fasta_target);
	my $tblastn_outfile;
	if(($output_fmt eq 'tab') or ($output_fmt eq 'all')){
		my $tblastn_outfile = $tblastn_filename . ".tsv.txt";
		unless(-s $tblastn_outfile){
			warn "Generating tblastn tab-delimited file....\n";
			my $tblastnCmd  = "$tblastn -query $fasta_query -db $fasta_target -max_target_seqs $num_alignments -evalue 10 -outfmt '6 qseqid salltitles qcovhsp pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
			warn $tblastnCmd . "\n\n";
			
			open(OUTFILE, ">$tblastn_outfile") or die "Couldn't open file $tblastn_outfile for writting, $!";
			print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "percent_positives", "align_length", "num_mismatch", 
			"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n"; 
			local (*TBLASTN_OUT, *TBLASTN_IN);
			my  $pid = open2(\*TBLASTN_OUT,\*TBLASTN_IN, $tblastnCmd) or die "Error calling open2: $!";
			close TBLASTN_IN or die "Error closing STDIN to tblastn process: $!";
			while(<TBLASTN_OUT>){
				chomp $_;
				my @blastn_hit =  split(/\t/, $_);
				my ($query_name, $target_name, $query_coverage, $percent_identity, $percent_positives, $align_length, $num_mismatch,
					$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @blastn_hit;
				if($percent_identity >= $min_percent_id){
					$e_value = "< 1e-179" if ($e_value =~ m/0\.0/);
					print OUTFILE join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $percent_positives, $align_length, $num_mismatch, 
						$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
				}
			}
			close TBLASTN_OUT or die "Error closing STDOUT from tblastn process: $!";
			wait;
			close(OUTFILE) or die "Couldn't close file $tblastn_outfile";
		}

	}
	
	if(($output_fmt eq 'align') or ($output_fmt eq 'all')){
		my $tblastn_outfile = $tblastn_filename . ".aln.txt";
		unless(-s $tblastn_outfile){
			warn "Generating tblastn alignment file....\n";
			my $tblastnCmd  = "$tblastn -query $fasta_query -db $fasta_target -num_descriptions $num_descriptions -num_alignments $num_alignments -evalue 10 -out $tblastn_outfile -num_threads $blast_num_cpu";
			warn $tblastnCmd . "\n\n";

			my $status = system($tblastn, 
				'-query', $fasta_query,
				'-db', $fasta_target,
				'-seg', 'yes',
				'-num_descriptions', $num_descriptions,
				'-num_alignments', $num_alignments,
				'-evalue', 10,
				'-out', $tblastn_outfile,
				'-num_threads', $blast_num_cpu
			) == 0 or die "Error calling $tblastn: $?";

		}
	}
	return $tblastn_outfile;
}
