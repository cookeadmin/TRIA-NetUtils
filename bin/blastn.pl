#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;

my ($target_infile, $query_infile, $min_percent_id, $blastn_task, $gap_open, $gap_extend, $num_descriptions, $num_alignments, $blast_num_cpu, $output_fmt, $output_dir);
GetOptions(
      'd=s'    => \$target_infile,
      'i=s'    => \$query_infile,
      'p=s'    => \$min_percent_id,
      't=s'    => \$blastn_task,
      'g=s'    => \$gap_open,
      'e=s'    => \$gap_extend,
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
$blastn_task = 'blastn' unless defined $blastn_task;
$gap_open = 5 unless defined $gap_open;
$gap_extend = 2 unless defined $gap_extend;
$num_descriptions = 25 unless defined $num_descriptions;
$num_alignments = 25 unless defined $num_alignments;
$blast_num_cpu = 2 unless defined $blast_num_cpu;
$output_fmt = 'all';

my ($makeblastdb, $blastn);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$blastn				= '/usr/local/bin/blastn';

sub usage {

die <<"USAGE";

Usage: $0 -d target_infile -i query_infile -p min_percent_id -t blastn_task -g gap_open -e gap_extend -v num_descriptions -b num_alignments -a blast_num_cpu -f output_fmt -o output_dir

Description - 

OPTIONS:

      -d target_infile -

      -i query_infile -

      -p min_percent_id 

      -t blastn_task - 

      -g gap_open - gap opening penalty (cost to open a gap).

      -e gap_extend - gap extension penalty (cost to extend a gap).

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
my $blastn_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . ".$blastn_task"));

my $blastn_infile = generate_blastn($query_infile, $target_infile, $blastn_task, $gap_open, $gap_extend, $num_descriptions, $num_alignments, $min_percent_id, $blast_num_cpu, $output_fmt, $blastn_filename);

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

sub generate_blastn{
	my $fasta_query = shift;
	die "Error lost fasta query file" unless defined $fasta_query;
	my $fasta_target = shift;
	die "Error lost fasta database target file" unless defined $fasta_target;
	my $blastn_task = shift;
	die "Error lost blastn program task" unless defined $blastn_task;
	my $gap_open = shift;
	die "Error lost gap opening penalty" unless defined $gap_open;
	my $gap_extend = shift;
	die "Error lost gap extension penalty" unless defined $gap_extend;
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
	my $blastn_filename = shift;
	die "Error lost blastn output filename" unless defined $blastn_filename;

	makeblastdb_nuc($fasta_target);

	my $blastn_outfile;
	if(($output_fmt eq 'tab') or ($output_fmt eq 'all')){
		my $blastn_outfile = $blastn_filename . ".tsv";
		unless(-s $blastn_outfile){
			warn "Generating blastn tab-delimited file....\n";
			my $blastnCmd  = "$blastn -query $fasta_query -db $fasta_target -task $blastn_task -gapopen $gap_open -gapextend $gap_extend -dust yes -max_target_seqs $num_alignments -evalue 1e-6 -outfmt '6 qseqid salltitles qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
			warn $blastnCmd . "\n\n";
			
			open(OUTFILE, ">$blastn_outfile") or die "Couldn't open file $blastn_outfile for writting, $!";
			print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
			"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n"; 
			local (*BLASTN_OUT, *BLASTN_IN);
			my  $pid = open2(\*BLASTN_OUT,\*BLASTN_IN, $blastnCmd) or die "Error calling open2: $!";
			close BLASTN_IN or die "Error closing STDIN to blastn process: $!";
			while(<BLASTN_OUT>){
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
			close BLASTN_OUT or die "Error closing STDOUT from blastn process: $!";
			wait;
			close(OUTFILE) or die "Couldn't close file $blastn_outfile";
		}

	}

	if(($output_fmt eq 'align') or ($output_fmt eq 'all')){
		my $blastn_outfile = $blastn_filename . ".aln";
		unless(-s $blastn_outfile){
			warn "Generating blastn alignment file....\n";
			my $blastnCmd  = "$blastn -query $fasta_query -db $fasta_target -task $blastn_task -gapopen $gap_open -gapextend $gap_extend -dust yes -num_descriptions $num_descriptions -num_alignments $num_alignments -evalue 1e-6 -out $blastn_outfile -num_threads $blast_num_cpu";
			warn $blastnCmd . "\n\n";

			my $status = system($blastn, 
				'-query', $fasta_query,
				'-db', $fasta_target,
				'-task', $blastn_task,
				'-gapopen', $gap_open,
				'-gapextend', $gap_extend,
				'-dust', 'yes',
				'-num_descriptions', $num_descriptions,
				'-num_alignments', $num_alignments,
				'-evalue', 1e-6,
				'-out', $blastn_outfile,
				'-num_threads', $blast_num_cpu
			) == 0 or die "Error calling $blastn: $?";

		}
	}
	return $blastn_outfile;
}
