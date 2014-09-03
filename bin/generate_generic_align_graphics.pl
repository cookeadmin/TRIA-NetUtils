#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use Bio::Graphics;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;

my ($fasta_infile, $blast_infile, $output_dir);
GetOptions(
      'i=s'    => \$fasta_infile,
      'b=s'    => \$blast_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $fasta_infile
      and defined $blast_infile
      and defined $output_dir
);

sub usage {

die <<"USAGE";

Usage: $0 -i fasta_infile -b blast_infile -o output_dir

Description - 

OPTIONS:

      -i fasta_infile -

      -b blast_infile -

      -o output_dir -

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %fasta_seq_lengths = ();
my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $sequence_length = length($sequence);
      $fasta_seq_lengths{$seq_id} = $sequence_length;
}


my %blast_entries = ();
open(INFILE, "<$blast_infile") or die "Couldn't open file $blast_infile for reading, $!";
my $i = 0;
while(<INFILE>){
      chomp $_;
      if($i ne 0){
#  	    warn $_ . "\n";
	    
	    my @split_row_entry = split(/\t/, $_);
	    my ($query_name,$target_name,$percent_identity,$align_length,$num_mismatch,
		  $num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score) = @split_row_entry;

	    my $blast_entry = join("\t",$target_name,$percent_identity,$align_length,$num_mismatch,
		  $num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score);

	    push(@{$blast_entries{$query_name}}, $blast_entry);
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $blast_infile";

foreach my $query_name (sort keys %blast_entries){

      warn "Processing blast results for " . $query_name . "....\n";

      my $query_length = $fasta_seq_lengths{$query_name};

      my $panel = Bio::Graphics::Panel->new(
	    -length    => $query_length,
	    -width     => 1000,
	    -pad_left  => 10,
	    -pad_right => 200,
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
	    -fgcolor     => 'black',
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

      foreach my $blast_entry (@{$blast_entries{$query_name}}){
	    my @split_blast_entry = split(/\t/, $blast_entry);
	    my ($target_name,$percent_identity,$align_length,$num_mismatch,
		  $num_gaps,$query_start,$query_end,$target_start,$target_end,$e_value,$bit_score) = @split_blast_entry;
	    my $strand_direction = "";
	    if($query_start < $query_end){
		  $strand_direction = "+";
	    }else{
		  $strand_direction = "-";
	    }
	    
	    if($query_start > $query_end){
		  my $temp_position = $query_start;
		  $query_start = $query_end;
		  $query_end = $temp_position;
	    }

	    my $strand = "";
	    if($strand_direction eq "+"){
		  $strand = 1;
	    }elsif($strand_direction eq "-"){
		  $strand = -1;
	    }

	    my ($target_id, $target_annotation) = "";
	    if($target_name =~ /\s\|\s/){
# 		  warn $target_name . "\n";
		  my @split_target_name = split(/\s\|\s/, $target_name);
		  $target_id = $split_target_name[0];
		  $target_annotation = $split_target_name[2];
	    }

	    my $description = join("; ", $target_annotation, join("=", "percent_idenity", $percent_identity), join("=", "length", $align_length), join(" = ", "strand", $strand_direction), join("=", "bit_score", $bit_score));
	    my $feature = Bio::SeqFeature::Generic->new(
		  -start        => $query_start, 
		  -end          => $query_end,
		  -score        => $bit_score,
		  -display_name => $target_id,
		  -tag          => {
			description => $description
		  },
	    );
	    
	    $track->add_feature($feature);
      }

      warn "Generating alignment image for " . $query_name . "....\n";
      my $outfile = join('/', $output_dir, join("-", $query_name, "alignment" . ".png"));
      open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";                          
      print OUTFILE $panel->png;
      close(OUTFILE) or die "Couldn't close file $outfile";
}

# my $searchio = Bio::SearchIO->new(-file   => $blast_infile,
#                                   -format => 'blast') or die "Error: Can't parse blast file $blast_infile using SearchIO $!";
# while(my $result = $searchio->next_result){
#       next unless($result->num_hits >= 1); # skip if there are no hits
# 
#       warn "Processing blast results for " . $result->query_name . "....\n";
# 
#       my $panel = Bio::Graphics::Panel->new(
# 	    -length    => $result->query_length,
# 	    -width     => 1000,
# 	    -pad_left  => 10,
# 	    -pad_right => 200,
#       );
#       
#       my $full_length = Bio::SeqFeature::Generic->new(
# 	    -start        => 1,
# 	    -end          => $result->query_length,
# 	    -display_name => $result->query_name,
#       );
# 
#       $panel->add_track(
# 	    $full_length,
# 	    -glyph   => 'arrow',
# 	    -tick    => 2,
# 	    -fgcolor => 'black',
# 	    -double  => 1,
# 	    -label   => 1,
#       );
#       
#       my $track = $panel->add_track(
# 	    -glyph       => 'graded_segments',
# 	    -label       => 1,
# 	    -connector   => 'dashed',
# 	    -bgcolor     => 'blue',
# 	    -font2color  => 'red',
# 	    -sort_order  => 'high_score',
# 	    -description => sub {
# 		  my $feature = shift;
# 		  return unless $feature->has_tag('description');
# 		  my $description = $feature->each_tag_value('description');
# 		  my $score = $feature->score;
# 		  "$description, score=$score";
# 	    },
#       );
#       
#       while( my $hit = $result->next_hit ) {
# 	    next unless(defined($hit));
# # 	    next unless $hit->significance < 1E-20;
# 	    my $feature = Bio::SeqFeature::Generic->new(
# 		  -score        => $hit->raw_score,
# 		  -display_name => $hit->name,
# 		  -tag          => {
# 			description => $hit->description
# 		  },
# 	    );
# 
# 	    while( my $hsp = $hit->next_hsp ) {
# 		  next unless(defined($hsp));
# 		  $feature->score($hsp->bits);
# 		  $feature->add_sub_SeqFeature($hsp,'EXPAND');
# 	    }
#       
# 	    $track->add_feature($feature);
#       }
#       
#       
#       my $outfile = join('/', $output_dir, join("-", $result->query_name, "alignment" . ".png"));
#       open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";                          
#       print OUTFILE $panel->png;
#       close(OUTFILE) or die "Couldn't close file $outfile";
# 
#       
# }

# Clean out the sequence I/O object.
$seqio = ();