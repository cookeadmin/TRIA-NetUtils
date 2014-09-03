#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;
use Bio::SeqIO;
use Switch;

# perl process_gmap_output.pl -i /home/cookeadmin/workspace/cathy_snps/AllSNPs_output_GMAP_2014-05-01.txt -o /home/cookeace/cathy_snps

my ($gmap_infile, $output_dir);
GetOptions(
      'i=s'    => \$gmap_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $gmap_infile
      and defined $output_dir
);

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i gmap_infile -o output_dir
    
Description - 
    
OPTIONS:
     
      -i gmap_infile - 
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %parse_gmap = ();
my $gmap_query_id = "";
open(INFILE, "<$gmap_infile") or die "Couldn't open file $gmap_infile for reading, $!";
while(<INFILE>){
      chomp $_;
# 	    warn "$_\n";
      $_ =~ s/^\s+//g;

      if($_ =~ m/^>([\w\_\-\.]+)$/){

	    $gmap_query_id = $1;
	    push(@{$parse_gmap{$gmap_query_id}}, $gmap_query_id);

      }elsif($_ !~ m/Alignments:/ and $_ !~ m/Alignment for path \d+:/ and $_ !~ m/^$/){
	    push(@{$parse_gmap{$gmap_query_id}}, $_);
      }
      
# 	    $gmap_contents .= $_ . "\n";
}
close(INFILE) or die "Couldn't close file $gmap_infile";

my $gmap_outfile = join('/', $output_dir, "gmap_outfile.txt");
open(OUTFILE, ">$gmap_outfile") or die "Couldn't open file $gmap_outfile for writting, $!";                          
print OUTFILE join("\t", "alignment_name","num_paths","coverage","percent_id","matches","mismatches","indels",
		  "unknowns","query_id","query_start","query_end","query_length","query_strand","target_id","target_start",
		  "target_end","target_length","target_strand","amino_acid_start","amino_acid_end","amino_acid_length",
		  "amino_acid_changes","num_exons","query_align_block","target_align_block", "intron_length_block") . "\n";
foreach my $gmap_query_id (keys %parse_gmap){
      my $gmap_entry = join("\t", @{$parse_gmap{$gmap_query_id}});

      my ($alignment_paths,$coverage,$percent_identity,$matches,$mismatches,$indels,
	    $unknowns,$query_id,$query_start,$query_end,$query_length,$query_strand,$target_id,
	    $target_start,$target_end,$target_length,$target_strand,$amino_acid_start,
	    $amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$intron_blocks);
      
# 	    my (@exon_align_identity,@query_align_blocks,@target_align_blocks,@intron_length_blocks) = ();
      my (@query_align_blocks,@target_align_blocks,@intron_length_blocks) = ();

      foreach my $field (@{$parse_gmap{$gmap_query_id}}){
	    
	    if($field =~ m/^([\w\_\-\.]+)$/){
		  $query_id = $1;
	    }elsif($field =~ m/Paths \((\d+)\):/){
		  
		  $alignment_paths = $1;
# 		  if($field =~ m/\*\*\*/){
# 			die "hello";
# 			$alignment_paths = "Possible chimera (Paths $1) ";
# 		  }
	    }elsif($field =~ m/Path \d+: query ([\d,]+)\.\.([\d,]+) \(([\d,]+) bp\) => genome ([\w\_\-\.]+):([\d,]+)\.\.([\d,]+) \(-*([\d,]+) bp\)/){
	    
		  ($query_start, $query_end, $query_length, $target_id, $target_start, $target_end, $target_length) = ($1, $2, $3, $4, $5, $6, $7);

		  $target_start =~ s/,//g;
		  $target_end =~ s/,//g;

		  my $target_temp = -1;
		  if($target_start > $target_end){
			$target_temp = $target_start;
			$target_start = $target_end;
			$target_end = $target_temp;
		  }

	    }elsif($field =~ m/cDNA direction: (sense|antisense|indeterminate)/){
		  
		  my $strandedness = $1; 
		  switch ($strandedness) {
			case "sense" {
			      $query_strand = "+"; 
			}
			case "antisense" { 
			      $query_strand = "-";
			}
			case "indeterminate" {
			      $query_strand = "?"; 
			}
			else { 
			      die "Strandedness: $strandedness" 
			}
		  }
	    }elsif($field =~ m/Genomic pos: [\w\_\-\.]+:[\d,]+,[\d,]+\.\.[\d,]+,[\d,]+ \(([+-]) strand\)/){
		  $target_strand = $1;
	    }elsif($field =~ m/Number of exons: (\d+)/){
		  $num_exons = $1;
	    }elsif($field =~ m/Coverage: (\d+\.\d+) \(query length: \d+ bp\)/){
		  $coverage = $1;
	    }elsif($field =~ m/Percent identity: (\d+\.\d+) \((\d+) matches, (\d+) mismatches, (\d+) indels, (\d+) unknowns\)/){

		  ($percent_identity,$matches,$mismatches,$indels,$unknowns) = ($1,$2,$3,$4,$5);

	    }elsif($field =~ m/Translation: (\d+)\.\.(\d+) \((\d+) aa\)/){

		  ($amino_acid_start,$amino_acid_end,$amino_acid_length) = ($1,$2,$3);

		  my $amino_acid_temp = -1;
		  if($amino_acid_start > $amino_acid_end){
			$amino_acid_temp = $amino_acid_start;
			$amino_acid_start = $amino_acid_end;
			$amino_acid_end = $amino_acid_temp;
		  }

	    }elsif($field =~ m/Amino acid changes: (([ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]) \[(\d+)\]|)/){
		  
		  my $has_aa_change = $2;
		  if(defined($has_aa_change)){
			my ($aa_change, $aa_position) = ($2,$3);
			$amino_acid_changes = join(";", $aa_change, $aa_position);
		  }else{
			$amino_acid_changes = "N/A";
		  }
		  
	    }
	    elsif($field =~ m/^[\+-][\w\_\-\.]+:(\d+-\d+)  \((\d+-\d+)\)   (\d+)%( (<-|\(-|==|-\)|->)   ...(\d+)...  (\d+\.\d+), (\d+\.\d+)|)/){
		  my ($target_align, $query_align, $align_identity, $has_intron) = ($1,$2,$3,$4);
# 			warn "$target_align, $query_align, $align_identity, $has_intron";
		  
# 			push(@exon_align_identity, $align_identity);

		  push(@target_align_blocks, $target_align);

		  push(@query_align_blocks, $query_align);

		  if(defined($has_intron) and $has_intron ne ""){
			      my ($intron_length, $intron_coverage, $intron_identity) = ($6,$7,$8);
			      $intron_blocks = join(";", $intron_length, $intron_coverage, $intron_identity);
# 			$intron_length = $6;


		  }else{

			$intron_blocks = "N/A";
		  }
		  push(@intron_length_blocks, $intron_blocks);
	    }
      }

# 	    my $exon_align = join(",", @exon_align_identity);
      my $query_align_block = join(",", @query_align_blocks);
      my $target_align_block = join(",", @target_align_blocks);
      my $intron_length_block = join(",", @intron_length_blocks);

# 	    $query_id =~ s/_/-/g;
# 	    $target_id =~ s/_/-/g;
	    print OUTFILE join("\t", join("-cmp-", $query_id, $target_id), $alignment_paths,$coverage,$percent_identity,$matches,$mismatches,$indels,
		  $unknowns,$query_id,$query_start,$query_end,$query_length,$query_strand,$target_id,
		  $target_start,$target_end,$target_length,$target_strand,$amino_acid_start,
		  $amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,
		  $query_align_block,$target_align_block, $intron_length_block) . "\n";

#       print OUTFILE join("\t", join("-vs-", $query_id, $target_id), $alignment_paths,$coverage,$percent_identity,$matches,$mismatches,$indels,
# 	    $unknowns,$query_id,$query_start,$query_end,$query_length,$query_strand,$target_id,
# 	    $target_start,$target_end,$target_length,$target_strand,$amino_acid_start,
# 	    $amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,
# 	    $query_align_block,$target_align_block, $intron_length_block) . "\n";
}
close(OUTFILE) or die "Couldn't close file $gmap_outfile";
