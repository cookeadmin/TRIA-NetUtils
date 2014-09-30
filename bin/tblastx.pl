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


$min_percent_id = 80 unless defined $min_percent_id;
$num_descriptions = 5 unless defined $num_descriptions;
$num_alignments = 5 unless defined $num_alignments;
$blast_num_cpu = 2 unless defined $blast_num_cpu;
$output_fmt = 'all';

my ($makeblastdb, $tblastx);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$tblastx			= '/usr/local/bin/tblastx';

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
my $tblastx_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . '.tblastx'));

my $tblastx_infile = generate_tblastx($query_infile, $target_infile, $num_descriptions, $num_alignments, $min_percent_id, $blast_num_cpu, $output_fmt, $tblastx_filename);

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
	my $output_fmt = shift;
	die "Error lost blastn output format" unless defined $output_fmt;
	my $tblastx_filename = shift;
	die "Error lost tblastx output filename" unless defined $tblastx_filename;

	makeblastdb_nuc($fasta_target);
	my $tblastx_outfile;
	if(($output_fmt eq 'tab') or ($output_fmt eq 'all')){
		my $tblastx_outfile = $tblastx_filename . ".tsv";
		unless(-s $tblastx_outfile){
			warn "Generating tblastx tab-delimited file....\n";
			my $tblastxCmd  = "$tblastx -query $fasta_query -db $fasta_target -seg yes -max_target_seqs $num_alignments -evalue 1e-6 -outfmt '6 qseqid salltitles qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
			warn $tblastxCmd . "\n\n";
			
			open(OUTFILE, ">$tblastx_outfile") or die "Couldn't open file $tblastx_outfile for writting, $!";
			print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
			"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n"; 
			local (*TBLASTX_OUT, *TBLASTX_IN);
			my  $pid = open2(\*TBLASTX_OUT,\*TBLASTX_IN, $tblastxCmd) or die "Error calling open2: $!";
			close TBLASTX_IN or die "Error closing STDIN to tblastx process: $!";
			while(<TBLASTX_OUT>){
				chomp $_;
				my @blastn_hit =  split(/\t/, $_);
				my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
					$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @blastn_hit;
				if($percent_identity >= $min_percent_id){
					$e_value = "< 1e-179" if ($e_value =~ m/0\.0/);
					print OUTFILE join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch, 
						$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
				}
			}
			close TBLASTX_OUT or die "Error closing STDOUT from tblastx process: $!";
			wait;
			close(OUTFILE) or die "Couldn't close file $tblastx_outfile";
		}

	}
	
	if(($output_fmt eq 'align') or ($output_fmt eq 'all')){
		my $tblastx_outfile = $tblastx_filename . ".aln";
		unless(-s $tblastx_outfile){
			warn "Generating tblastx alignment file....\n";
			my $tblastxCmd  = "$tblastx -query $fasta_query -db $fasta_target -seg yes -num_descriptions $num_descriptions -num_alignments $num_alignments -evalue 1e-6 -out $tblastx_outfile -num_threads $blast_num_cpu";
			warn $tblastxCmd . "\n\n";

			my $status = system($tblastx, 
				'-query', $fasta_query,
				'-db', $fasta_target,
				'-seg', 'yes',
				'-num_descriptions', $num_descriptions,
				'-num_alignments', $num_alignments,
				'-evalue', 1e-6,
				'-out', $tblastx_outfile,
				'-num_threads', $blast_num_cpu
			) == 0 or die "Error calling $tblastx: $?";

		}
	}
	return $tblastx_outfile;
}
