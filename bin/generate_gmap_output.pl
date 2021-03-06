#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;

# A program for automating gmap_build and gmap to align cDNA sequences to a set of Genome sequences.
my ($db_name, $gmap_target_file, $gmap_query_file, $kmer_size, $gmap_num_cpu, $gmap_num_ram, $output_dir);
GetOptions(
      'n=s'    => \$db_name, # The name of the Genome database
      't=s'    => \$gmap_target_file, # The target genome input file
      'q=s'    => \$gmap_query_file, # The query cdna input file
      'k=s'    => \$kmer_size, # The k-mer value for genomic index
      'c=s'    => \$gmap_num_cpu, # The number of cpu threads to allocate to gmap
      'm=s'    => \$gmap_num_ram, # The amount of ram to allocate to gmap
      'o=s'    => \$output_dir, # The output directory for gmap results
);

usage() unless (
      defined $db_name
      and defined $gmap_target_file
      and defined $gmap_query_file
      and defined $output_dir
);

$kmer_size = 15 unless defined $kmer_size; # k-mer size is set to 15 as default unless specified by input parameters
$gmap_num_cpu = 2 unless defined $gmap_num_cpu; # number of cpus is set to 2 as default unless specified by input parameters
$gmap_num_ram = 3 unless defined $gmap_num_ram; # amount of ram is set to 3 as default unless specified by input parameters

# Paths to the gmap_build and gmap program executables.
my ($gmap_build, $gmap);
$gmap_build			= '/usr/local/bin/gmap_build';
$gmap				= '/usr/local/bin/gmap';

sub usage {

die <<"USAGE";

Usage: $0 -n db_name -t gmap_target_file -q gmap_query_file -k kmer_size -c gmap_num_cpu -m gmap_num_ram -o output_dir

Description - 

OPTIONS:
      -n db_name - The name of the Genome database

      -t gmap_target_file - The target genome input file

      -q gmap_query_file - The query cdna input file

      -k kmer_size - The k-mer value for genomic index

      -c gmap_num_cpu - The number of cpu threads to allocate to gmap

      -m gmap_num_ram - The amount of ram to allocate to gmap

      -o output_dir - The output directory for gmap results

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Generate gmap output alignment summary output file.
generate_gmap($db_name, $gmap_target_file, $gmap_query_file, $kmer_size, $output_dir, $gmap_num_cpu, $gmap_num_ram);


#gmap_build -D ./ -d Pinus_taeda pita_v1.01_maker.20131025.151125.fsa -k 15
sub gmap_build{

      my $db_name = shift;
      die "Error lost db_name in sub gmap_build" unless defined $db_name;

      my $gmap_target_file = shift;
      die "Error lost gmap target genome input file in sub gmap_build" unless defined $gmap_target_file;
      
      my $kmer_size = shift;
      die "Error lost k-mer value for genomic index in sub gmap_build" unless defined $kmer_size;

      my $output_dir = shift;
      die "Error lost output_dir in sub gmap_build" unless defined $output_dir;
      
      

      my $db_dir = join('/', $output_dir, $db_name);

      


	    warn "Calling gmap_build for $gmap_target_file....\n";
	    # gmap_build -D ./ -d Pinus_taeda pita_v1.01_maker.20131025.151125.fsa -k 15
	    warn "$gmap_build -D $output_dir -d $db_name $gmap_target_file -k $kmer_size\n";
	    # die join(", ", $db_name, $gmap_target_file, $kmer_size, $output_dir);
	    system($gmap_build, 
		  '-D', $output_dir, 
		  '-d', $db_name,
		  $gmap_target_file,
		  '-k', $kmer_size
	    ) == 0 or die "Error calling $gmap_build -D $output_dir -d $db_name $gmap_target_file -k $kmer_size: $?";
      return $db_dir;
}

# gmap -D ./Pinus_taeda -d Pinus_taeda -S -A -5 AllSNPsContigs.fa -t 2 -B 3 > gmap_output_test.txt
sub generate_gmap{
 
      my $db_name = shift;
      die "Error lost db_name in sub generate_gmap" unless defined $db_name;

      my $gmap_target_file = shift;
      die "Error lost gmap target genome input file in sub generate_gmap" unless defined $gmap_target_file;

      my $gmap_query_file = shift;
      die "Error lost gmap query cdna input file in sub generate_gmap" unless defined $gmap_query_file;

      my $kmer_size = shift;
      die "Error lost k-mer value for genomic index in sub generate_gmap" unless defined $kmer_size;

      my $output_dir = shift;
      die "Error lost output_dir in sub generate_gmap" unless defined $output_dir;

      my $gmap_num_cpu = shift;
      die "Error lost the number of cpu threads in sub generate_gmap" unless defined $gmap_num_cpu;

      my $gmap_num_ram = shift;
      die "Error lost the number of ram in sub generate_gmap" unless defined $gmap_num_ram;


      my $gmap_outfile = join("_", $gmap_target_file, $gmap_query_file . '.gmap');

      my $db_dir = gmap_build($db_name, $gmap_target_file, $kmer_size, $output_dir);

      warn "Generating gmap output file....\n";
       my $gmapCmd  = "$gmap -D $db_dir -d $db_name -S -A $gmap_query_file -O -t $gmap_num_cpu -B $gmap_num_ram";
#      my $gmapCmd  = "$gmap -D $db_dir -d $db_name -S $gmap_query_file -n 1 -O -t $gmap_num_cpu -B $gmap_num_ram";

      warn $gmapCmd . "\n\n";

      local (*GMAP_OUT, *GMAP_IN);
      my  $pid = open2(\*GMAP_OUT,\*GMAP_IN, $gmapCmd) or die "Error calling open2 for $gmap process: $!";
      close GMAP_IN or die "Error closing STDIN to $gmap process blastn: $!";	

      open(OUTFILE, ">$gmap_outfile") or die "Couldn't open file $gmap_outfile for writting, $!";
      #print OUTFILE join("\t", "query_name", "target_name", "percent_identity", "align_length", "num_mismatch", 
	#    "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
      while (<GMAP_OUT>){ 
	    chomp $_;
	    print OUTFILE $_ . "\n";
      }
      close(OUTFILE) or die "Couldn't close file $gmap_outfile";

      close GMAP_OUT or die "Error closing STDOUT from $gmap process: $!";
      wait;

      return $gmap_outfile;
}
