#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
use Switch;

# perl process_gmap_output-new.pl -i /home/cookeadmin/workspace/cathy/AllSNPs_output_GMAP -o /home/cookeadmin/cathy_snps

my ($gmap_infile, $query_infile, $output_dir);
GetOptions(
      'i=s'    => \$gmap_infile,
      'q=s'    => \$query_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $gmap_infile
      and defined $query_infile
      and defined $output_dir
);

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i gmap_infile -q query_infile -o output_dir
    
Description - 
    
OPTIONS:
     
      -i gmap_infile - 
      -q query_infile - 
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %query_unique_ids;
my $seqio = Bio::SeqIO->new(-file => $query_infile, '-format' => 'Fasta');
while(my $seq_entry = $seqio->next_seq) {

      my $seq_id = $seq_entry->id;
      my $sequence = $seq_entry->seq;
      my $seq_desc = $seq_entry->desc;

      my $fasta_unique_id;
      if($seq_desc eq ""){
	    $fasta_unique_id = $seq_id;
      }else{
	    
	   $fasta_unique_id = join(" ", $seq_id, $seq_desc);
      }
      $query_unique_ids{$fasta_unique_id} = $fasta_unique_id;
}

my %parsed_gmap = ();
my $gmap_query_id = "";
open(INFILE, "<$gmap_infile") or die "Couldn't open file $gmap_infile for reading, $!";
while(<INFILE>){
      chomp $_;
# 	    warn "$_\n";
      $_ =~ s/^\s+//g;

      if($_ =~ m/^>(\w|\d)/){
	    $gmap_query_id = $_;
	    $gmap_query_id =~ s/>//g;
	    if(defined($query_unique_ids{$gmap_query_id})){
# 		  warn $gmap_query_id . "\n";
		  push(@{$parsed_gmap{$gmap_query_id}}, $gmap_query_id);
	    }
      }elsif($_ !~ m/^\s+/ and $_ !~ m/^$/){

	    push(@{$parsed_gmap{$gmap_query_id}}, $_);
      }
      
}
close(INFILE) or die "Couldn't close file $gmap_infile";

my $outfilename = fileparse($gmap_infile);
my $gmap_outfile = join('/', $output_dir, $outfilename . ".parsed");
open(OUTFILE, ">$gmap_outfile") or die "Couldn't open file $gmap_outfile for writting, $!";                          
print OUTFILE join("\t", "alignment_name","coverage","percent_id","matches","mismatches","indels",
		  "unknowns","query_id","query_start","query_end","query_length","query_strand","target_id","target_start",
		  "target_end","target_length","target_strand","amino_acid_start","amino_acid_end","amino_acid_length",
		  "amino_acid_changes","num_exons","query_align_block","target_align_block", "align_identity_block", 
		  "intron_length_block","target_gap_positions","query_gap_positions","mismatch_positions") . "\n";
my %gmap_output_final = ();
foreach my $gmap_query_id (keys %parsed_gmap){
      my %processed_gmap = ();
      my $gmap_entry = join("\t", @{$parsed_gmap{$gmap_query_id}});

      my ($alignment_paths,$path_num,$coverage,$percent_identity,$matches,$mismatches,$indels,
	    $unknowns,$query_id,$query_start,$query_end,$query_length,$query_strand,$target_id,
	    $target_start,$target_end,$target_length,$target_strand,$amino_acid_start,
	    $amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$intron_blocks);

      my ($current_query_start, $current_query_end, $current_query_length, $current_target_id, 
      $current_target_start, $current_target_end, $current_target_length);
      my $current_path = 1;
      
      foreach my $field (@{$parsed_gmap{$gmap_query_id}}){
	    
	    if($field =~ m/^($gmap_query_id)/){
		  $query_id = $1;
	    }

	    if($field =~ m/Paths \((\d+)\):/){
		  
		  $alignment_paths = $1;
		  if($alignment_paths eq 0){ # Exit out of loop if $alignment_paths equals 0 indicating no alignment info to parse.
			last;
		  }
	    }
	    
	    if($field =~ m/Path (\d+): query ([\d,]+)\.\.([\d,]+) \(([\d,]+) bp\) => genome ([\w\_\-\.\d]+):([\d,]+)\.\.([\d,]+) \(-*([\d,]+) bp\)/){
	    
		  ($path_num, $query_start, $query_end, $query_length, $target_id, $target_start, $target_end, $target_length) = ($1, $2, $3, $4, $5, $6, $7, $8);

		  $target_start =~ s/,//g;
		  $target_end =~ s/,//g;

		  my $target_temp = -1;
		  if($target_start > $target_end){
			$target_temp = $target_start;
			$target_start = $target_end;
			$target_end = $target_temp;
		  }
	    }

	    
	    if(defined($path_num)){
		  warn "path: $path_num eq $current_path";
		  if($path_num eq $current_path){
			if($field =~ m/cDNA direction: (sense|antisense|indeterminate)/){
			      
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
			}
			if($field =~ m/Genomic pos: [\w\_\-\.\d]+:[\d,]+,[\d,]+\.\.[\d,]+,[\d,]+ \(([+-]) strand\)/){
			      $target_strand = $1;
			}
			if($field =~ m/Number of exons: (\d+)/){
			      $num_exons = $1;
			}
			if($field =~ m/Coverage: (\d+\.\d+) \(query length: \d+ bp\)/){
			      $coverage = $1;
			}
			if($field =~ m/Percent identity: (\d+\.\d+) \((\d+) matches, (\d+) mismatches, (\d+) indels, (\d+) unknowns\)/){

			      ($percent_identity,$matches,$mismatches,$indels,$unknowns) = ($1,$2,$3,$4,$5);

			}
			if($field =~ m/Translation: (\d+)\.\.(\d+) \((\d+) aa\)/){

			      ($amino_acid_start,$amino_acid_end,$amino_acid_length) = ($1,$2,$3);

			      my $amino_acid_temp = -1;
			      if($amino_acid_start > $amino_acid_end){
				    $amino_acid_temp = $amino_acid_start;
				    $amino_acid_start = $amino_acid_end;
				    $amino_acid_end = $amino_acid_temp;
			      }

			}
			if($field =~ m/Amino acid changes: ((.+)|)/){
			      
			      my $has_aa_change = $2;
			      if(defined($has_aa_change)){
				    my $amino_acid_changes_list = $2;
				    my @amino_acid_change_entries = split(/, /, $amino_acid_changes_list);
				    my @amino_acid_changes_final = ();
				    foreach my $amino_acid_change_entry (@amino_acid_change_entries){
					  warn $amino_acid_change_entry . "\n";
					  my ($aa_change, $aa_position);
					  if($amino_acid_change_entry =~ m/(.+) \[(\d+)\]/){
						($aa_change, $aa_position) = ($1, $2);
					  }
					  push(@amino_acid_changes_final, join(";", $aa_change, $aa_position));
				    }

				    $amino_acid_changes = join(",", @amino_acid_changes_final);

			      }else{
				    $amino_acid_changes = "N/A";
			      }
			      ($current_query_start, $current_query_end, $current_query_length, $current_target_id, $current_target_start, $current_target_end, $current_target_length) = ($query_start, $query_end, $query_length, $target_id, $target_start, $target_end, $target_length);
			}
		  }
		  if($path_num > $current_path){

      # 		  print OUTFILE join("\t", join("-", $query_id, "cmp", $current_target_id, join("_", "path", $current_path, "of", $alignment_paths)),$coverage,$percent_identity,$matches,$mismatches,$indels,
      # 		  $unknowns,$query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
      # 		  $target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons) . "\n";
			$amino_acid_start = "N/A" unless(defined($amino_acid_start));
			$amino_acid_end = "N/A" unless(defined($amino_acid_end));
			$amino_acid_length = "N/A" unless(defined($amino_acid_length));

			$processed_gmap{$current_path} = join("\t", join("-", $query_id, "cmp", $current_target_id, "path" . $current_path . "of" . $alignment_paths),$coverage,$percent_identity,$matches,$mismatches,$indels,
			$unknowns,$query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
			$target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons);
			$current_path++;

			warn "path: $path_num > $current_path";
		  }
	    }

	    
      }
      next if ($alignment_paths eq 0); # Exit out of loop if $alignment_paths equals 0 indicating no alignment info to parse.
#       print OUTFILE join("\t", join("-", $query_id, "cmp", $current_target_id, join("_", "path", $current_path, "of", $alignment_paths)),$coverage,$percent_identity,$matches,$mismatches,$indels,
#       $unknowns,$query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
#       $target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons) . "\n";

      $amino_acid_start = "N/A" unless(defined($amino_acid_start));
      $amino_acid_end = "N/A" unless(defined($amino_acid_end));
      $amino_acid_length = "N/A" unless(defined($amino_acid_length));

      $processed_gmap{$current_path} = join("\t", join("-", $query_id, "cmp", $current_target_id, "path" . $current_path . "of" . $alignment_paths),$coverage,$percent_identity,$matches,$mismatches,$indels,
      $unknowns,$query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
      $target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons);
			
      my $found_align_tag = "false";
      my $found_target_align = "false";
      my (@query_align_blocks,@target_align_blocks,@intron_length_blocks,@align_identity_blocks,@query_alignment,@target_alignment,@consensus_alignment) = ();
      
      my ($block_path_num,$query_align_block,$target_align_block,$intron_length_block,$align_identity_block,$query_align_sequence,$target_align_sequence,$consensus_align_sequence);

      my $current_block_path = 1;
      foreach my $field (@{$parsed_gmap{$gmap_query_id}}){
	    if($field =~ m/Alignments:/){
		  $found_align_tag = "true";
# 		  die $found_align_tag;
	    }

	    if($found_align_tag eq "true"){
		  
		  if($field =~ m/Alignment for path (\d+):/){
			$block_path_num = $1;
		  }

		  if(defined($block_path_num)){
			warn "Block path: $block_path_num eq $current_block_path";
			if($block_path_num eq $current_block_path){
			      if($field =~ m/^[\+-][\w\_\-\.\d]+:(\d+-\d+)  \((\d+-\d+)\)   (\d+)%( (<-|\(-|==|-\)|->)   ...(\d+)...  (\d+\.\d+), (\d+\.\d+)|)/){
				    my ($target_align, $query_align, $align_identity, $has_intron) = ($1,$2,$3,$4);
		  # 			warn "$target_align, $query_align, $align_identity, $has_intron";

				    push(@target_align_blocks, $target_align);

				    push(@query_align_blocks, $query_align);
				    
				    push(@align_identity_blocks, $align_identity);
				    if(defined($has_intron) and $has_intron ne ""){
					  my ($intron_length, $intron_coverage, $intron_identity) = ($6,$7,$8);
					  $intron_blocks = join(";", $intron_length, $intron_coverage, $intron_identity);
				    }else{

					  $intron_blocks = "N/A";
				    }
				    push(@intron_length_blocks, $intron_blocks);
			      }

			      #  +tscaffold7014:134446 TAGCTTGACAATGTAAAAAAGCCTCAGAATTTCATGCTTTTGTCTATTAC
			      if($field =~ m/^[\+-][\w\_\-\.\d]+:\d+\s([AGCTRYKMSWBDHVN\.\s]+)/i  and $found_target_align eq "false" and ($field !~ m/^\d+\s{2,}$/ and $field !~ m/^\d+\s{5}/)){
				    my $target_align_partial = $1;
				    push(@target_alignment, $target_align_partial);
				    $found_target_align = "true";
			      }

			      # 46 TAGCTTGACAATGTAAAAAAGCCTCAGAATTTCATGCTTTTGTCTATTAC
			      if($field =~ m/^\d+\s([AGCTRYKMSWBDHVN\s\d]+)/i and $found_target_align eq "true"){
				    print $field . "\n";
				    my $query_align_partial = $1;
				    push(@query_alignment, $query_align_partial);
				    $found_target_align = "false";
			      }

			}
			if($block_path_num > $current_block_path){

			      $query_align_block = join(",", @query_align_blocks);
			      $target_align_block = join(",", @target_align_blocks);
			      $intron_length_block = join(",", @intron_length_blocks);
			      $align_identity_block = join(",", @align_identity_blocks);
			      
			      $query_align_sequence = join("", @query_alignment);
			      $target_align_sequence = join("", @target_alignment);
			      
			      $target_align_sequence =~ s/[AGCTRYKMSWBDHVN]{3}\.\.\.[AGCTRYKMSWBDHVN]{3}//gi;
			      $target_align_sequence =~ s/\s/-/g;
			      $query_align_sequence =~ s/\s{4}\d{1}\s{4}//g;
			      $query_align_sequence =~ s/\s{3}\d{2}\s{4}//g;
			      $query_align_sequence =~ s/\s{3}\d{3}\s{3}//g;
			      $query_align_sequence =~ s/\s{2}\d{4}\s{3}//g;
			      $query_align_sequence =~ s/\s{2}\d{5}\s{2}//g;
			      $query_align_sequence =~ s/\s{1}\d{6}\s{2}//g;
			      $query_align_sequence =~ s/\s{1}\d{7}\s{1}//g;
			      $query_align_sequence =~ s/\s/-/g;

			      my ($query_align_seq_length,$target_align_seq_length) = (length($query_align_sequence),length($target_align_sequence));
			      warn join("\t",$processed_gmap{$current_block_path});
			      print $target_align_sequence . "\n";
			      print $query_align_sequence . "\n";
			      die "Error: query_align_seq_length=$query_align_seq_length and target_align_seq_length=$target_align_seq_length are not equal" unless($target_align_seq_length eq $query_align_seq_length);

			      warn join("\t", $target_align_seq_length, $query_align_seq_length) . "\n";

			      my @target_align_seq_chars = split('', $target_align_sequence);
			      my @query_align_seq_chars = split('', $query_align_sequence);

			      my (@target_gap_positions,@query_gap_positions,@mismatch_positions) = ();
			      for(my $i = 0; $i < $query_align_seq_length; $i++){
				    my $position = ($i + 1);
				    if($target_align_seq_chars[$i] ne $query_align_seq_chars[$i]){
					  warn join("\t", $target_align_seq_chars[$i], $query_align_seq_chars[$i]) . "\n";
					  if($target_align_seq_chars[$i] eq "-"){
						push(@target_gap_positions, $position);
					  }elsif($query_align_seq_chars[$i] eq "-"){
						push(@query_gap_positions, $position);
					  }else{
						push(@mismatch_positions, $position);
					  }
				    }
			      }

			      my ($target_gap_position_list,$query_gap_position_list,$mismatch_position_list) = "";
			      $target_gap_position_list = join(",", @target_gap_positions);
			      $query_gap_position_list = join(",", @query_gap_positions);
			      $mismatch_position_list = join(",", @mismatch_positions);

			      $target_gap_position_list = "N/A" if($target_gap_position_list eq "");
			      $query_gap_position_list = "N/A" if($query_gap_position_list eq "");
			      $mismatch_position_list = "N/A" if($mismatch_position_list eq "");

			      $processed_gmap{$current_block_path} = join("\t",$processed_gmap{$current_block_path},$query_align_block,$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,$query_gap_position_list,$mismatch_position_list);
			      (@query_align_blocks,@target_align_blocks,@intron_length_blocks,@align_identity_blocks,@query_alignment,@target_alignment,@consensus_alignment) = ();

			      warn $processed_gmap{$current_block_path} . "\n";
			      warn $target_align_sequence . "\n";
			      warn $query_align_sequence . "\n";

			      warn "Block path: $block_path_num > $current_block_path";
			      $current_block_path++;
			}
		  }
	    }
      }

      $query_align_block = join(",", @query_align_blocks);
      $target_align_block = join(",", @target_align_blocks);
      $intron_length_block = join(",", @intron_length_blocks);
      $align_identity_block = join(",", @align_identity_blocks);

      $query_align_sequence = join("", @query_alignment);
      $target_align_sequence = join("", @target_alignment);

      $target_align_sequence =~ s/[AGCTRYKMSWBDHVN]{3}\.\.\.[AGCTRYKMSWBDHVN]{3}//gi;
      $target_align_sequence =~ s/\s/-/g;
      $query_align_sequence =~ s/\s{4}\d{1}\s{4}//g;
      $query_align_sequence =~ s/\s{3}\d{2}\s{4}//g;
      $query_align_sequence =~ s/\s{3}\d{3}\s{3}//g;
      $query_align_sequence =~ s/\s{2}\d{4}\s{3}//g;
      $query_align_sequence =~ s/\s{2}\d{5}\s{2}//g;
      $query_align_sequence =~ s/\s{1}\d{6}\s{2}//g;
      $query_align_sequence =~ s/\s{1}\d{7}\s{1}//g;
      $query_align_sequence =~ s/\s/-/g;


      my ($query_align_seq_length,$target_align_seq_length) = (length($query_align_sequence),length($target_align_sequence));
      warn join("\t",$processed_gmap{$current_block_path});
      print $target_align_sequence . "\n";
      print $query_align_sequence . "\n";

      die "Error: query_align_seq_length=$query_align_seq_length and target_align_seq_length=$target_align_seq_length are not equal" unless($target_align_seq_length eq $query_align_seq_length);

      warn join("\t", $target_align_seq_length, $query_align_seq_length) . "\n";

      my @target_align_seq_chars = split('', $target_align_sequence);
      my @query_align_seq_chars = split('', $query_align_sequence);

      my (@target_gap_positions,@query_gap_positions,@mismatch_positions) = ();
      for(my $i = 0; $i < $query_align_seq_length; $i++){
	    my $position = ($i + 1);
	    if($target_align_seq_chars[$i] ne $query_align_seq_chars[$i]){
		  warn join("\t", $target_align_seq_chars[$i], $query_align_seq_chars[$i]) . "\n";
		  if($target_align_seq_chars[$i] eq "-"){
			push(@target_gap_positions, $position);
		  }elsif($query_align_seq_chars[$i] eq "-"){
			push(@query_gap_positions, $position);
		  }else{
			push(@mismatch_positions, $position);
		  }
	    }
      }

      my ($target_gap_position_list,$query_gap_position_list,$mismatch_position_list) = "";
      $target_gap_position_list = join(",", @target_gap_positions);
      $query_gap_position_list = join(",", @query_gap_positions);
      $mismatch_position_list = join(",", @mismatch_positions);

      $target_gap_position_list = "N/A" if($target_gap_position_list eq "");
      $query_gap_position_list = "N/A" if($query_gap_position_list eq "");
      $mismatch_position_list = "N/A" if($mismatch_position_list eq "");

      $processed_gmap{$current_block_path} = join("\t",$processed_gmap{$current_block_path},$query_align_block,$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,$query_gap_position_list,$mismatch_position_list);
      warn $processed_gmap{$current_block_path} . "\n";
      warn $target_align_sequence . "\n";
      warn $query_align_sequence . "\n";

      foreach my $processed_path_num (keys %processed_gmap){
	    my @processed_gmap_entry = split(/\t/, $processed_gmap{$processed_path_num});
	    my $query_name = $processed_gmap_entry[7];
	    $gmap_output_final{$query_name}{$processed_path_num} = $processed_gmap{$processed_path_num};
      }


}

foreach my $query_name (sort {$a cmp $b} keys %gmap_output_final){
      foreach my $processed_path_num (sort %{$gmap_output_final{$query_name}}){
	    print OUTFILE $gmap_output_final{$query_name}{$processed_path_num} . "\n" if(defined($gmap_output_final{$query_name}{$processed_path_num}));
      }
}
close(OUTFILE) or die "Couldn't close file $gmap_outfile";
