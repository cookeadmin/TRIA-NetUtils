#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

use Getopt::Long;

# perl annotate_heatmap_R_input.pl -l /home/cookeadmin/workspace/adriana/Phloem/phloem-microarray-annotations/blastx_output/phloem-blastx-file-list -m /home/cookeadmin/workspace/adriana/Phloem/heatmap_input_R-2014-06-18/JP-DE-all-days-comparison/submission-summary/JP-DE-all-days-comparison.tsv -f /home/cookeadmin/workspace/adriana/ncbi_nr_flat_files_2014-07-02/ncbi_nr_db_notes.txt -o /home/cookeadmin/workspace/adriana/Phloem/annotated_heatmap_input_R-2014-06-27/JP-DE-all-days-comparison

my ($mapman_id_list, $gene_exp_matrix_infile, $output_dir);
my @options = (
      'l=s'    => \$mapman_id_list,
      'm=s'    => \$gene_exp_matrix_infile,
      'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
      defined $blastx_annotation_list
      and $gene_exp_matrix_infile
      and $output_dir
);

sub usage {
    
    die << "USAGE";
    
Usage: $0 -i mapman_id_list -m gene_exp_matrix_infile -o output_dir
    
USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %gene_exp_matrix_id_list = ();
open(INFILE, "<$gene_exp_matrix_infile") or die "Couldn't open file $gene_exp_matrix_infile for reading, $!";
my $i = 0;
while(<INFILE>){
      chomp $_;
#       warn $_ . "\n";
      if($i eq 0){
	    $gene_exp_matrix_id_list{"HEADER"} = $_;

      }else{
	    my @split_row_entry = split(/\t/, $_);
	    my $annotation = $split_row_entry[0];
	    my @split_annotation_entry = split(/ /, $annotation, 2);
	    my $name_id = $split_annotation_entry[0];
	    $gene_exp_matrix_id_list{$name_id} = $_;
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $gene_exp_matrix_infile";

my %mapman_entry_list = ();
open(INFILE, "<$mapman_id_list") or die "Couldn't open file $mapman_id_list for reading, $!";
while(<INFILE>){
      chomp $_;
       warn $_ . "\n";
      my $name_id = $_;
      $mapman_entry_list{$name_id} = $name_id;
}
close(INFILE) or die "Couldn't close file $mapman_id_list";

my $gene_exp_matrix_filename = fileparse($gene_exp_matrix_infile);
my $mapman_id_list_filename  = fileparse($mapman_id_list);
my $outfile = join('/', $output_dir, join("_", $mapman_id_list_filename, $gene_exp_matrix_filename));
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE $gene_exp_matrix_id_list{"HEADER"} . "\n";
foreach my $name_id (sort {$a cmp $b} keys %blast_annotation_list){
      print OUTFILE $gene_exp_matrix_id_list{$name_id} . "\n" if(defined($gene_exp_matrix_id_list{$name_id}));
}
close(OUTFILE) or die "Couldn't close file $outfile";