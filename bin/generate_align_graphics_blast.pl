#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::Graphics;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;
# use Data::Dump qw(dump);

# perl generate_align_graphics.pl -i ~/workspace/loblloy_pine_microarray_analysis/jack_pine_blastx/jack_pine_illumina_trinity_23Aug2011.fasta_jack-pine-illumina-trinity-23Aug2011-megablast-loblloly-microarray-contigs-list.fasta_TAIR10_pep_20101214_updated.blastx.aln -o ~/workspace/loblloy_pine_microarray_analysis/jack_pine_blastx/alignment_images_blast
my ($blast_infile, $output_dir);
GetOptions(
      'i=s'    => \$blast_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $blast_infile
      and defined $output_dir
);

sub usage {

die <<"USAGE";

Usage: $0 -i blast_infile -o output_dir

Description - 

OPTIONS:

      -i blast_infile - The default blast output file as input. (-outfmt 0 = pairwise).

      -o output_dir - The output directory to write the Bio::Graphics blast alignment images.

USAGE
}


# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Array to store blast hit information for each query.
my @blast_results = ();

warn "Processing blast results for " . $blast_infile . "....\n";
my $searchio = Bio::SearchIO->new(-file   => $blast_infile,
                                  -format => 'blast') or die "Error: Can't parse blast file $blast_infile using SearchIO $!";
while(my $result = $searchio->next_result){
      next unless($result->num_hits >= 1);
      my $query_name = $result->query_name;
      my $query_length = $result->query_length;
      my $blast_algorithm = $result->algorithm;

      warn "Generating alignment image for " . $query_name . "....\n";
      while(my $hit = $result->next_hit){
	    while(my $hsp = $hit->next_hsp){

		  my ($target_name,$target_header,$percent_identity,$align_length,$num_mismatch,
			$num_gaps,$blast_frame,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score);

		  $target_name = $hit->name;

		  $target_header = $hit->description;

		  $percent_identity = sprintf("%.2f", $hsp->percent_identity);

		  $align_length= $hsp->length('total');

		  $num_mismatch = ($hsp->length('total') - ($hsp->num_identical + $hsp->gaps));

		  $num_gaps = $hsp->gaps;

		  $blast_frame = (($hsp->query->frame + 1) * $hsp->query->strand);
		  if($blast_frame >= 0){
			$query_start = $hsp->start('query');
			$query_end = $hsp->end('query');
		  }elsif($blast_frame < 0){
			$query_start = $hsp->end('query');
			$query_end = $hsp->start('query');
		  }

		  $target_start = $hsp->start('hit');
		  $target_end = $hsp->end('hit');

		  $e_value = $hsp->evalue;
		  $e_value = "< 1e-179" if ($e_value =~ m/0\.0/);

		  $bit_score = $hsp->bits;

		  my $blast_results_entry = join("\t", $target_name,$target_header,$percent_identity,$align_length,
			$num_mismatch, $num_gaps,$blast_frame,$query_start,$query_end,$target_start,$target_end,
			$e_value,$bit_score);

		  warn "$query_name ==> $target_name" . "\n";
		  push(@blast_results, $blast_results_entry);
		  
	    }
      }

      my $panel = Bio::Graphics::Panel->new(
	    -length    => $query_length,
	    -width     => 1000,
	    -pad_left  => 10,
	    -pad_right => 600,
      );
      
      my $full_length = Bio::SeqFeature::Generic->new(
	    -start        => 1,
	    -end          => $query_length,
	    -display_name => $query_name,
      );

      $panel->add_track(
	    $full_length,
	    -glyph   => 'arrow',
	    -tick    => 2,
	    -fgcolor => 'black',
	    -double  => 1,
	    -label   => 1,
      );
      
      my $track = $panel->add_track(
	    -glyph       => 'graded_segments',
	    -label       => 1,
	    -connector   => 'dashed',
	    -bgcolor     => 'blue',
	    -font2color  => 'red',
	    -sort_order  => 'high_score',
	    -description => sub {
			my $feature = shift;
			return unless $feature->has_tag('description');
			my $description = $feature->each_tag_value('description');
      # 		  my $score = $feature->score;
			"$description";
	    },
      );

      foreach my $blast_result_entry (@blast_results) {

	    my @split_blast_result_entry = split(/\t/, $blast_result_entry);
	    my ($target_name,$target_header,$percent_identity,$align_length,$num_mismatch,$num_gaps,
		  $blast_frame,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score) = @split_blast_result_entry;

	    if($query_start > $query_end){
		  my $temp_position = $query_start;
		  $query_start = $query_end;
		  $query_end = $temp_position;
	    }

	    if($blast_frame > 0){
		  $blast_frame = "+" . $blast_frame;
	    }

	    my $target_annotation = "";
	    if($target_header =~ /\s\|\s/){
# 		  warn $target_name . "\n";
		  my @split_target_header = split(/\s\|\s/, $target_header);
		  $target_annotation = $split_target_header[1];
	    }

	    if($blast_algorithm eq "BLASTN"){
		  $align_length = join(" ", $align_length, "nt");
	    }elsif($blast_algorithm eq "BLASTX"){
		  $align_length = join(" ", $align_length, "aa");
	    }

	    my $description = join("; ", $target_annotation, join("= ", "percent_idenity", $percent_identity), join("= ", "align_length", $align_length), join("= ", "frame", $blast_frame), join("= ", "e_value", $e_value), join("= ", "bit_score", $bit_score));
	    warn join(", ", $query_name, $target_name, $description) . "\n";
	    my $feature = Bio::SeqFeature::Generic->new(
		  -start        => $query_start, 
		  -end          => $query_end,
		  -score        => $bit_score,
		  -display_name => $target_name,
		  -tag          => {
			description => $description
		  },
	    );
	    
	    $track->add_feature($feature);
      }    

      
      my $outfile = join('/', $output_dir, join("-", $query_name, "alignment" . ".png"));
      open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";                          
      print OUTFILE $panel->png;
      close(OUTFILE) or die "Couldn't close file $outfile";

      @blast_results = (); # Clear the array after every query entry
}

