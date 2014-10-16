#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;

my ($target_infile, $query_infile, $blast_num_cpu, $output_fmt, $output_dir);
GetOptions(
      'd=s'    => \$target_infile,
      'i=s'    => \$query_infile,
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
$output_fmt = 'all' unless defined $output_fmt;

my ($makeblastdb, $blastn);

$makeblastdb 			= '/usr/bin/makeblastdb';
$blastn				= '/usr/bin/blastn';

sub usage {

die <<"USAGE";

Usage: $0 -d target_infile -i query_infile -a blast_num_cpu -f output_fmt -o output_dir

Description - 

OPTIONS:

      -d target_infile -

      -i query_infile -

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
my $primer_blastn_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . ".primer_blastn"));

my $primer_blastn_infile = generate_primer_blastn($query_infile, $target_infile, $blast_num_cpu, $output_fmt, $primer_blastn_filename);

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

sub generate_primer_blastn{
	my $fasta_query = shift;
	die "Error lost fasta query file" unless defined $fasta_query;
	my $fasta_target = shift;
	die "Error lost fasta database target file" unless defined $fasta_target;
	my $blast_num_cpu = shift;
	die "Error lost number of cpus to allocate" unless defined $blast_num_cpu;
	my $output_fmt = shift;
	die "Error lost primer blastn output format" unless defined $output_fmt;
	my $primer_blastn_outfile = shift;
	die "Error lost primer blastn output filename" unless defined $primer_blastn_outfile;

	makeblastdb_nuc($fasta_target);

	my $primer_blastn_outfile;
	if(($output_fmt eq 'tab') or ($output_fmt eq 'all')){
		my $primer_blastn_outfile = $primer_blastn_filename . ".tsv";
		unless(-s $primer_blastn_outfile){
			warn "Generating blastn tab-delimited file....\n";
			my $primerBlastnCmd  = "$blastn -query $fasta_query -db $fasta_target -task blastn -word_size 7 -dust no -evalue 1000 -outfmt '6 qseqid salltitles qcovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
			warn $primerBlastnCmd . "\n\n";
			
			open(OUTFILE, ">$primer_blastn_outfile") or die "Couldn't open file $primer_blastn_outfile for writting, $!";
			print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", 
			"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n"; 
			local (*PRIMER_BLASTN_OUT, *PRIMER_BLASTN_IN);
			my $pid = open2(\*PRIMER_BLASTN_OUT,\*PRIMER_BLASTN_IN, $primerBlastnCmd) or die "Error calling open2: $!";
			close PRIMER_BLASTN_IN or die "Error closing STDIN to primer blastn process: $!";
			while(<PRIMER_BLASTN_OUT>){
				chomp $_;
				my @primer_blastn_hit =  split(/\t/, $_);
				my ($query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch,
					$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @blastn_hit;
				if($percent_identity >= $min_percent_id){
					$e_value = "< 1e-179" if ($e_value =~ m/0\.0/);
					print OUTFILE join("\t", $query_name, $target_name, $query_coverage, $percent_identity, $align_length, $num_mismatch, 
						$num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
				}
			}
			close PRIMER_BLASTN_OUT or die "Error closing STDOUT from primer blastn process: $!";
			wait;
			close(OUTFILE) or die "Couldn't close file $primer_blastn_outfile";
		}

	}

	if(($output_fmt eq 'align') or ($output_fmt eq 'all')){
		my $primer_blastn_outfile = $primer_blastn_filename . ".aln";
		unless(-s $primer_blastn_outfile){
			warn "Generating primer blastn alignment file....\n";
			my $primerBlastnCmd  = "$blastn -query $fasta_query -db $fasta_target -task blastn -word_size 7 -dust no -evalue 1000 -out $primer_blastn_outfile -num_threads $blast_num_cpu";
			warn $primerBlastnCmd . "\n\n";

			my $status = system($blastn, 
				'-query', $fasta_query,
				'-db', $fasta_target,
				'-task', 'blastn',
				'-word_size', 7,
				'-dust', 'no',
				'-evalue', 1000,
				'-out', $primer_blastn_outfile,
				'-num_threads', $blast_num_cpu
			) == 0 or die "Error calling $blastn: $?";

		}
	}
	return $primer_blastn_outfile;
}
