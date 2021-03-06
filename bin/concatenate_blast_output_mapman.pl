#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

use Getopt::Long;

# perl concatenate_blast_output_mapman.pl -l /home/cookeadmin/workspace/adriana/Phloem/phloem-microarray-annotations/blastx_output-2014-07-03/phloem-blastx-file-list -m /home/cookeadmin/workspace/adriana/Phloem/Phloem-comparisons-annotated.csv -f /home/cookeadmin/workspace/adriana/ncbi_nr_flat_files_2014-07-02/ncbi_nr_db_notes.txt -o /home/cookeadmin/workspace/test
# perl concatenate_blast_output_mapman.pl -l /home/cookeadmin/workspace/adriana/Xylem/xylem-microarray-annotations/blastx_output-2014-07-03/xylem-blastx-file-list -m /home/cookeadmin/workspace/adriana/Xylem/ -f /home/cookeadmin/workspace/adriana/ncbi_nr_flat_files_2014-07-02/ncbi_nr_db_notes.txt -o /home/cookeadmin/workspace/test
my ($blastx_annotation_list, $gene_exp_matrix_infile, $ncbi_flat_notes_infile, $output_dir);
my @options = (
      'l=s'    => \$blastx_annotation_list,
      'm=s'    => \$gene_exp_matrix_infile,
      'f=s'    => \$ncbi_flat_notes_infile,
      'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
      defined $blastx_annotation_list
      and $gene_exp_matrix_infile
      and $ncbi_flat_notes_infile
      and $output_dir
);

sub usage {
    
    die << "USAGE";
    
Usage: $0 -l blastx_annotation_list -m gene_exp_matrix_infile -f ncbi_flat_notes_infile -o output_dir
    
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
#        warn $_ . "\n";
      if($i eq 0){
	    $gene_exp_matrix_id_list{"HEADER"} = $_;

      }else{
	    my @split_row_entry = split(/\t/, $_);
	    my $name_id = $split_row_entry[6];
	    push(@{$gene_exp_matrix_id_list{$name_id}}, $_);
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $gene_exp_matrix_infile";

my %ncbi_flat_notes_list = ();
open(INFILE, "<$ncbi_flat_notes_infile") or die "Couldn't open file $ncbi_flat_notes_infile for reading, $!";
$i = 0;
while(<INFILE>){
      chomp $_;
#        warn $_ . "\n";
      if($i ne 0){
	    my @split_row_entry = split(/\t/, $_);
	    my ($target_name, $target_notes) = ($split_row_entry[0], $split_row_entry[1]);

	    my @split_target_notes = split("\001", $target_notes);
	    my $target_notes_desc = $split_target_notes[0];
	    $ncbi_flat_notes_list{$target_name} = $target_notes_desc;
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $ncbi_flat_notes_infile";

my @blast_annotation_list = ();
open(INFILE, "<$blastx_annotation_list") or die "Couldn't open file $blastx_annotation_list for reading, $!";
while(<INFILE>){
      chomp $_;
#       warn $_ . "\n";
      my @split_row_entry = split(/\t/, $_);
      my ($blast_db, $blast_input_dir, @blast_files) = ($split_row_entry[0], $split_row_entry[1], @split_row_entry[2..$#split_row_entry]);
      
      foreach my $blast_file (@blast_files){
	    my $blast_filename = join('/', $blast_input_dir, $blast_file);
	    push(@blast_annotation_list, $blast_filename);
	    
      }

}
close(INFILE) or die "Couldn't close file $blastx_annotation_list";

my %blast_annotations = ();
my %name_id_list = ();

foreach my $blast_filename (@blast_annotation_list){

      warn $blast_filename . "\n";
      open(INFILE, "<$blast_filename") or die "Couldn't open file $blast_filename for reading, $!";
      my $i = 0;
      while(<INFILE>){
	    chomp $_;
	    if($i ne 0){
# 			warn $_ . "\n";
		  
		  my @split_row_entry = split(/\t/, $_);
		  my $query_id = $split_row_entry[0];
		  my @split_query_id = split("\001", $query_id);
		  my $name_id = $split_query_id[0];

		  # fix < 1e-179 so that we can sort by evalue numerically
		  if($split_row_entry[10] =~ m/< 1e-179/){
			$split_row_entry[10] = "1e-179";
		  }

		  my $blast_entry = join("\t", @split_row_entry);
# 			warn $blast_entry . "\n";

		  push(@{$blast_annotations{$name_id}}, [split(/\t/, $blast_entry)]);
	    }
	    $i++;
      }
      close(INFILE) or die "Couldn't close file $blast_filename";
}
    


my %blast_final_list = ();
foreach my $name_id (sort keys %blast_annotations){
# 	    warn $name_id . "\n";
      if(defined(@{$blast_annotations{$name_id}})){
      
	    my @blast_annotation_sorted = sort {$b->[11] <=> $a->[11]} @{$blast_annotations{$name_id}};
      
	    foreach my $blast_annotation_entry (@blast_annotation_sorted){
		  $blast_final_list{$name_id} = join("\t", @$blast_annotation_entry);
		  warn join("\t", @$blast_annotation_entry) . "\n";
 		  last;
	    }

      }else{
	    next;
      }
}

my %gene_exp_matrix_final = ();
foreach my $name_id (sort keys %gene_exp_matrix_id_list){
# 	    warn $name_id . "\n";

      if($name_id ne "HEADER"){
	    foreach my $gene_exp_matrix_id_list_entry (@{$gene_exp_matrix_id_list{$name_id}}){
			
			my $annotated_matrix_id_entry = "";
			if(defined($blast_final_list{$name_id})){
			      my @split_blast_final_list = split(/\t/, $blast_final_list{$name_id});
			      my ($query_name,$target_name,$percent_identity,$align_length,$num_mismatch,$num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score) = @split_blast_final_list;

			      my $target_annotation = "";
			      if($target_name =~ m/\001/){
				    my @split_target_name = split("\001", $target_name);
				    $target_annotation = $split_target_name[0];
			      }elsif($target_name !~ m/gi\|\d+\|/){
				    my @split_target_name = split(/ \| /, $target_name);
				    my ($target_id, $target_desc) = ($split_target_name[0], $split_target_name[2]);
				    if($target_desc =~ /;/){
					  my @split_target_desc = split(/;/, $target_desc);
					  $target_desc = $split_target_desc[0];
				    }
				    if($target_desc =~ /\//){
					  my @split_target_desc = split(/\//, $target_desc);
					  $target_desc = $split_target_desc[0];
				    }

				    $target_annotation = join(" ", $target_id, $target_desc, "[Arabidopsis thaliana]");

			      }else{
				    $target_annotation = $target_name;
			      }
			      
			      if($target_annotation =~ /unknown/ and defined($ncbi_flat_notes_list{$target_name})){
				    my $new_annotation = $ncbi_flat_notes_list{$target_name};
				    $new_annotation = substr($new_annotation, 0, 50) . "..." if(length($new_annotation) >= 50);
				    $target_annotation =~ s/unknown/$new_annotation/;
			      }
		  # 		  die join("\t", $query_name,$target_name,$percent_identity,$align_length,$num_mismatch,$num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score);
			      push(@{$gene_exp_matrix_final{$name_id}}, join("\t", $gene_exp_matrix_id_list_entry, $query_name, $target_annotation, $bit_score));

			}
		  
	    }
      }
}



my $filename = fileparse($gene_exp_matrix_infile);
my $outfile = join('/', $output_dir, $filename);
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE join("\t", $gene_exp_matrix_id_list{"HEADER"}, "query_id", "target_id", "bit_score") . "\n";
foreach my $name_id (keys %gene_exp_matrix_final){
      foreach my $gene_exp_matrix_id_list_entry (@{$gene_exp_matrix_final{$name_id}}){
	    print OUTFILE $gene_exp_matrix_id_list_entry . "\n";
      }
}
close(OUTFILE) or die "Couldn't close file $outfile";