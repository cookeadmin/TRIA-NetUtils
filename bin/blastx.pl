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

$blast_num_cpu = 2 unless defined $blast_num_cpu;
$min_percent_id = 80 unless defined $min_percent_id;
$num_descriptions = 5 unless defined $num_descriptions;
$num_alignments = 5 unless defined $num_alignments;
$output_fmt = 'all';

my ($makeblastdb, $blastx);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$blastx				= '/usr/local/bin/blastx';

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
my $blastx_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . '.blastx'));

my $blastx_infile = generate_blastx($query_infile, $target_infile, $min_percent_id, $num_descriptions, $num_alignments, $blastx_filename);

# makeblastdb -in ncbi_nr_db_2014-05-30_organisms.fasta -dbtype 'prot' -out ncbi_nr_db_2014-05-30_organisms.fasta
sub makeblastdb_pep{
      my $fastadb = shift;
      die "Error lost fastadb to makeblastdb" unless defined $fastadb;
      # format the database file into .pin .psq .phr files.
      my ($fastadbPIN, $fastadbPSQ, $fastadbPHR);
      $fastadbPIN = $fastadb . '.pin';
      $fastadbPSQ = $fastadb . '.psq';
      $fastadbPHR = $fastadb . '.phr';
      unless(-s $fastadbPIN and -s $fastadbPSQ and -s $fastadbPHR){
	    warn "Calling makeblastdb for $fastadb....\n";
	    warn "$makeblastdb -in $fastadb -dbtype prot\n\n";
	    system($makeblastdb, 
		  '-in', $fastadb, 
		  '-dbtype', 'prot'
	    ) == 0 or die "Error calling $makeblastdb -in $fastadb -dbtype prot: $?";
      }

}

sub generate_blastx{
	my $fasta_query = shift;
	die "Error lost fasta query file" unless defined $fasta_query;
	my $fasta_target = shift;
	die "Error lost fasta database target file" unless defined $fasta_target;
	my $min_percent_id = shift;
	die "Error lost minimum percent identity" unless defined $min_percent_id;
	my $num_descriptions = shift;
	die "Error lost number of descriptions" unless defined $num_descriptions;
	my $num_alignments = shift;
	die "Error lost number of alignments" unless defined $num_alignments;
	my $output_fmt = shift;
	die "Error lost blastn output format" unless defined $output_fmt;
	my $blastx_filename = shift;
	die "Error lost blastx output filename" unless defined $blastx_filename;

	makeblastdb_pep($fasta_target);
	
	my $blastx_outfile;
	if(($output_fmt eq 'tab') or ($output_fmt eq 'all')){
		my $blastx_outfile = $blastx_filename . ".tsv";
		unless(-s $blastx_outfile){
			warn "Generating blastx tab-delimited file....\n";
			my $blastxCmd  = "$blastx -query $fasta_query -db $fasta_target -seg no -max_target_seqs $num_alignments -evalue 1e-6 -outfmt '6 qseqid salltitles qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
			warn $blastxCmd . "\n\n";
			
			open(OUTFILE, ">$blastx_outfile") or die "Couldn't open file $blastx_outfile for writting, $!";
			print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
			"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n"; 
			local (*BLASTX_OUT, *BLASTX_IN);
			my  $pid = open2(\*BLASTX_OUT,\*BLASTX_IN, $tblastxCmd) or die "Error calling open2: $!";
			close BLASTX_IN or die "Error closing STDIN to blastx process: $!";
			while(<BLASTX_OUT>){
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
			close BLASTX_OUT or die "Error closing STDOUT from blastx process: $!";
			wait;
			close(OUTFILE) or die "Couldn't close file $blastx_outfile";
		}
	}
	if(($output_fmt eq 'align') or ($output_fmt eq 'all')){
		my $blastx_outfile = $blastx_filename . ".aln";
		unless(-s $blastx_outfile){
			warn "Generating blastx alignment file....\n";
			my $blastxCmd  = "$blastx -query $fasta_query -db $fasta_target -seg no -num_descriptions $num_descriptions -num_alignments $num_alignments -evalue 1e-6 -out $blastx_outfile -num_threads $blast_num_cpu";
			warn $blastxCmd . "\n\n";

			my $status = system($blastx, 
				'-query', $fasta_query,
				'-db', $fasta_target,
				'-seg', 'no',
				'-num_descriptions', $num_descriptions,
				'-num_alignments', $num_alignments,
				'-evalue', 1e-6,
				'-out', $blastx_outfile,
				'-num_threads', $blast_num_cpu
			) == 0 or die "Error calling $blastx: $?";


		}
	}
	return $blastx_outfile;
}
