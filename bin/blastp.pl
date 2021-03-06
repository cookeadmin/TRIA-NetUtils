#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use IPC::Open2;
# remember to fix join portion of script in tblastx, blastx, tblastn, and so on.
# perl blastp.pl -i /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GWAS/pgi/pgi_nr_protein_sequences/ncbi_nr_db-2014-09-10.fasta_insecta_taxa_list.txt.fasta_pgi_nr_proteins.txt.fasta -d /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GWAS/all_mpb_protein_sequences.fasta -p 0 -a 7 -v 1000000 -b 1000000 -o /home/cookeadmin/workspace/GBS_data-08-10-2013/MPB_GBS_Data-08-10-2013/MPB_GWAS/pgi/pgi_blastp

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


$output_fmt = 'all' unless defined $output_fmt;
$blast_num_cpu = 2 unless defined $blast_num_cpu;
$min_percent_id = 80 unless defined $min_percent_id;
$num_descriptions = 5 unless defined $num_descriptions;
$num_alignments = 5 unless defined $num_alignments;

my ($makeblastdb, $blastp);
$makeblastdb 			= '/usr/local/bin/makeblastdb';
$blastp				= '/usr/local/bin/blastp';

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
my $blastp_filename = join('/', $output_dir, join("_", $fasta_query_name, $fasta_target_name . '.blastp'));

my $blastp_infile = generate_blastp($query_infile, $target_infile, $min_percent_id, $num_descriptions, $num_alignments, $output_fmt, $blastp_filename);


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

sub generate_blastp{
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
	die "Error lost blastp output format" unless defined $output_fmt;
	my $blastp_filename = shift;
	die "Error lost blastp output filename" unless defined $blastp_filename;
    
	makeblastdb_pep($fasta_target);
	
	my $blastp_outfile;
	if(($output_fmt eq 'tab') or ($output_fmt eq 'all')){
		my $blastp_outfile = $blastp_filename . ".tsv.txt";
		unless(-s $blastp_outfile){
			warn "Generating blastp tab-delimited file....\n";
			my $blastpCmd  = "$blastp -query $fasta_query -db $fasta_target -seg no -max_target_seqs $num_alignments -evalue 1e-6 -outfmt '6 qseqid salltitles qcovhsp pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads $blast_num_cpu";
			warn $blastpCmd . "\n\n";
			
			open(OUTFILE, ">$blastp_outfile") or die "Couldn't open file $blastp_outfile for writting, $!";
			print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "percent_positives", "align_length", "num_mismatch",
			"num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
			local (*BLASTP_OUT, *BLASTP_IN);
			my  $pid = open2(\*BLASTP_OUT,\*BLASTP_IN, $blastpCmd) or die "Error calling open2: $!";
			close BLASTP_IN or die "Error closing STDIN to blastp process: $!";
			while(<BLASTP_OUT>){
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
			close BLASTP_OUT or die "Error closing STDOUT from blastp process: $!";
			wait;
			close(OUTFILE) or die "Couldn't close file $blastp_outfile";
		}
	}
	if(($output_fmt eq 'align') or ($output_fmt eq 'all')){
		my $blastp_outfile = $blastp_filename . ".aln.txt";
		unless(-s $blastp_outfile){
			warn "Generating blastp alignment file....\n";
			my $blastpCmd  = "$blastp -query $fasta_query -db $fasta_target -seg no -num_descriptions $num_descriptions -num_alignments $num_alignments -evalue 1e-6 -out $blastp_outfile -num_threads $blast_num_cpu";
			warn $blastpCmd . "\n\n";
            
			my $status = system($blastp,
				'-query', $fasta_query,
				'-db', $fasta_target,
				'-seg', 'no',
				'-num_descriptions', $num_descriptions,
				'-num_alignments', $num_alignments,
				'-evalue', 1e-6,
				'-out', $blastp_outfile,
				'-num_threads', $blast_num_cpu
			) == 0 or die "Error calling $blastp: $?";
            
            
		}
	}
	return $blastp_outfile;
}
