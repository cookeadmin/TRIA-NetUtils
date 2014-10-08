#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

use Getopt::Long;

# perl get_best_blast_hits.pl -l /home/cookeadmin/workspace/adriana/Phloem/phloem-microarray-annotations/blastx_output-2014-07-03/phloem-blastx-file-list -f /home/cookeadmin/workspace/adriana/ncbi_nr_flat_files_2014-07-02/ncbi_nr_db_notes.txt -n phloem-master-blastx-file -o /home/cookeadmin/workspace/test

my ($blastn_infile, $num_hits, $output_dir);
my @options = (
      'i=s'    => \$blastn_infile,
      'n=s'    => \$num_hits,
      'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
      defined $blastn_infile
      and $output_dir
);

$num_hits = 3 unless(defined($num_hits));

sub usage {
    
    die << "USAGE";
    
Usage: $0 -i blastn_infile -o output_dir
   
USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $blastn_filename = fileparse($blastn_infile, qr/\.[^.]*/);
my %blast_annotations = ();
open(INFILE, "<$blastn_infile") or die "Couldn't open file $blastn_infile for reading, $!";
my $i = 0;
while(<INFILE>){
      chomp $_;
      if($i ne 0){
# 	warn $_ . "\n";
	    
	    my @split_row_entry = split(/\t/, $_);
	    my $query_id = $split_row_entry[0];
	    my @split_query_id = split("\001", $query_id);
	    my $name_id = $split_query_id[0];
	    

	    my $blast_entry = join("\t", @split_row_entry);
# 			warn $blast_entry . "\n";

	    push(@{$blast_annotations{$name_id}}, [split(/\t/, $blast_entry)]);
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $blastn_infile";

   


my %blast_final_list = ();
foreach my $name_id (sort keys %blast_annotations){
# 	    warn $name_id . "\n";
      
      my $hits_counter = 1;
      if(defined(@{$blast_annotations{$name_id}})){
      
	    my @blast_annotation_sorted = sort {$b->[12] <=> $a->[12]} @{$blast_annotations{$name_id}};
      
	    foreach my $blast_annotation_entry (@blast_annotation_sorted){
		  
		  warn join("\t", @$blast_annotation_entry) . "\n";
		  if($hits_counter <= $num_hits){

			warn "$hits_counter <= $num_hits\n";
			push(@{$blast_final_list{$name_id}}, join("\t", @$blast_annotation_entry));
		  }
		  $hits_counter++;
	    }
      }else{
	    next;
      }
}

my %blast_entries = ();
foreach my $name_id (keys %blast_final_list){
      foreach my $blast_entry (@{$blast_final_list{$name_id}}){
	    # warn $name_id . "\n";
	    push(@{$blast_entries{$name_id}}, [split(/\t/, $blast_entry)]);
      }
}

# Get filename based on number of hits choosen.
my $outfile = "";
if($num_hits > 1){
	$outfile = join('/', $output_dir, join("_", $blastn_filename, "top", $num_hits, "hits") . ".txt");
}else{
	$outfile = join('/', $output_dir, join("_", $blastn_filename, "top", "hits") . ".txt");
}

open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "percent_identity", "align_length", "num_mismatch", "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
foreach my $name_id (sort {$a cmp $b} keys %blast_entries){
      foreach my $blast_entry (sort {$b->[12] <=> $a->[12]} @{$blast_entries{$name_id}}){
	    warn join("\t", @$blast_entry) . "\n";

	    my ($query_name,$target_name,$query_coverage,$percent_identity,$align_length,$num_mismatch,$num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score) = @$blast_entry;
	    print OUTFILE join("\t", $query_name,$target_name,$query_coverage,$percent_identity,$align_length,$num_mismatch,$num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score) . "\n";

      }
}
close(OUTFILE) or die "Couldn't close file $outfile";
