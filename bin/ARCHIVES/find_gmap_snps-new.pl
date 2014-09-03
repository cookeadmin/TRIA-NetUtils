#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;
use Bio::SeqIO;
use File::Basename;
use Switch;
# perl find_gmap_snps-new.pl -s /home/cookeadmin/workspace/cathy/AllSNPsVariantPositions.tsv -g /home/cookeadmin/workspace/cathy/AllSNPs_output_GMAP.parsed -q /home/cookeadmin/workspace/cathy/AllSNPsContigs.fa -o ~/workspace/cathy
my ($snps_infile, $gmap_infile, $query_infile, $output_dir);
GetOptions(
      's=s'    => \$snps_infile,
      'g=s'    => \$gmap_infile,
      'q=s'    => \$query_infile,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $snps_infile
      and defined $gmap_infile
      and defined $query_infile
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -s snps_infile -g gmap_infile -q query_infile -o output_dir
    
Description - 
    
OPTIONS:
      -s snps_infile - 
    
      -g gmap_infile - 

      -q query_infile -

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# One thing I forgot to mention, is that some of the SNPs may represent paralogues, this is confirmed when the SNP site lines up to more than one loblolly contig with high probability (one of the reasons for using the loblolly data to do this)
# 
# 1. Find SNP in GMAP output (if it is possible, record location of SNP if it is different from the position indicated by the SNP name)
# 
# 2. Does the SNP align to only one contig
#     NO: mark as potential paralogue go to next SNP
#     YES: Go to step 2
# 
# 2. Is the SNP in a coding region
#     NO: mark as 5' or 3' untranslated go to next SNP
#                 On a rare occasion, the SNP may be in the intron region, if that is the case, record as intronic and go to next SNP
#     YES: Go to step 3
# 
# 3. Identify the 3bp codon and the amino acid based on the direction of the cDNA (i.e. sense or anitsense)
# 
# 4. Look up the amino acid table - does the SNP variant result in an amino acid change?
#     NO: mark as synonymous mutation
#     YES: mark as non-synonymous and record the amino acid change


my $query_seqs = parse_fasta_infile($query_infile);

my $gmap_alignments = parse_gmap_infile($gmap_infile);

my $snps_positions = parse_snps_infile($snps_infile);


my (%paralogs, %coding_region, %five_prime_UTR, %three_prime_UTR, %indeterminate) = ();
foreach my $alignment_name (keys %{$gmap_alignments}){
      warn $alignment_name . "\n";
      my @split_gmap_align_entry = split(/\t/, $gmap_alignments->{$alignment_name});
      my ($coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,$query_start,
	    $query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,$target_strand,
	    $amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$query_align_block,
	    $target_align_block,$align_identity_block,$intron_length_block) = @split_gmap_align_entry;


#       warn join("\t", $alignment_name,$coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,$query_start,
# 	    $query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,$target_strand,
# 	    $amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$query_align_block,
# 	    $target_align_block,$align_identity_block,$intron_length_block) . "\n";

      #### When anti-sense works make sure that if ($num_align_paths > 1) and ($coverage > 90) and there are more than one 
      # found name them paralogs and exclude from this list if two or more paths meet the criteria and put into paralog file. 
      # if only one > 90 coverage then do not exclude and treat as non-paralogous snp.
      my ($path_num, $num_align_paths);
      if($alignment_name =~ m/path(\d+)of(\d+)/){
	    ($path_num, $num_align_paths) = ($1, $2);
# 	    warn join("\t", $path_num, $num_align_paths) . "\n";
      }

      
#       if(($num_align_paths > 1) and ($coverage > 90) and defined(@{$snps_positions->{$query_id}})){
# 	    warn $gmap_alignments->{$alignment_name} . "\n";
# 	    foreach my $snps_entry (@{$snps_positions->{$query_id}}){
# 		  warn $snps_entry . "\n";
# 		  my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);
# 
# 		  my $snp_base = get_subseq($query_seqs->{$query_id}, $snp_position, $snp_position);
# 		  my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
# 		  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
# 			$paralogs{$alignment_name} = join("\t",$alignment_name,$query_id,"SNP_Paralog",$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]"));
# 		  }else{
# 			$indeterminate{$alignment_name} = join("\t",$alignment_name,$query_id,"SNP_Paralog",$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]"));
# 		  }
# 		  
# 	    }
# 			      
#       ######   REVAMP PARALOGS LATER 
#       ##### When anti-sense strand works make >= 1 or leave out that clause so that we get the ones that arent excluded from 
#       # the paralog list. Should still take coverage > 90 so that we get good results.
#       }elsif(($num_align_paths > 1) and defined(@{$snps_positions->{$query_id}})){
#       if(($num_align_paths >= 1) and ($coverage >= 80) and defined(@{$snps_positions->{$query_id}})){
      if(($num_align_paths >= 1) and defined(@{$snps_positions->{$query_id}})){

	    if(($amino_acid_start ne "N/A") and ($amino_acid_end ne "N/A") and ($amino_acid_length ne "N/A")){
		  switch ($query_strand) {
			case "+" {
			      warn "IN SENSE STRAND (+)\n";
# 			      warn $gmap_alignments->{$alignment_name} . "\n";
			      warn $query_id . "\n";
			      my @query_subseq = split('', get_subseq($query_seqs->{$query_id}, $amino_acid_start, $amino_acid_end));

			      my $query_subsequence = join('', @query_subseq);

			      my $aa_seq = translate_dna($query_subsequence);
# 			      warn "$query_id\n" . join('', @query_subseq) . "\n";
# 			      warn "$query_id\n$aa_seq\n";
# 			      warn join("\t", $amino_acid_start, $amino_acid_end) . "\n";
			      foreach my $snps_entry (@{$snps_positions->{$query_id}}){
				    warn $snps_entry . "\n";
				    my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);

				    warn "($snp_position >= $query_start) and ($snp_position <= $query_end)" unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
				    next unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
				    if(($snp_position >= $amino_acid_start) and ($snp_position <= $amino_acid_end)){
					  my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence) = get_aa_codon_snp($query_subsequence,$amino_acid_start,$snp_position,$snp_variants);
					  
					  my @snp_variants_list_entries = split(/\//, $snp_variants_list);
					  my $num_snp_variants = scalar(@snp_variants_list_entries);
					  if($num_snp_variants eq 3){
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"CDS_SNP",$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence"));
					  }else{

						$coding_region{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"CDS_SNP",$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence"));
					  }
					  
				    }elsif($snp_position < $amino_acid_start){
					  warn "5 PRIME UTR\n";
					  my $five_prime_UTR_start = $target_start;
					  my $five_prime_UTR_end = $target_start + ($amino_acid_start - 1);
# 						my $five_prime_UTR_start = $target_start;
# 						my $five_prime_UTR_end = $target_end;

					  my $snp_base = get_subseq($query_seqs->{$query_id}, $snp_position, $snp_position);
					  my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
					  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
						$five_prime_UTR{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"five_prime_UTR_SNP",$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]"));
					  }elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
						my $snp_variant3 = $snp_base;
						my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"five_prime_UTR_SNP",$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]"));
					  }
				    }
				    elsif($snp_position > $amino_acid_end){
					  warn "3 PRIME UTR\n";
					  my $three_prime_UTR_start = $target_start + ($amino_acid_end + 1);
					  my $three_prime_UTR_end = $target_end;
# 						my $three_prime_UTR_start = $target_start;
# 						my $three_prime_UTR_end = $target_end;

					  my $snp_base = get_subseq($query_seqs->{$query_id}, $snp_position, $snp_position);
					  my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
					  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
						$three_prime_UTR{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"three_prime_UTR_SNP",($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]"));
					  }elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
						my $snp_variant3 = $snp_base;
						my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"three_prime_UTR_SNP",($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]"));
					  }
				    }
			      }
			}
			case "-" { 
			      warn "IN ANTISENSE STRAND (-)\n";
# 			      warn $gmap_alignments->{$alignment_name} . "\n";
			      warn $query_id . "\n";
# 			      warn $query_seqs->{$query_id} . "\n";
			      
			      my @query_subseq = split('', get_subseq($query_seqs->{$query_id}, $amino_acid_start, $amino_acid_end));
	    
	    
			      my $query_subsequence = join('', @query_subseq);
			      my $comp_query_subseq = reverse_complement(join('', @query_subseq));
	    # 
	    # 			# Since the cDNA direction is antisense the protein should not be reversed.
			      my $aa_seq = translate_dna($comp_query_subseq);
	    
# 			      warn "$query_id\n" . join('', @query_subseq) . "\n";
# 			      warn "$query_id\n$aa_seq\n";
			      
# 			      warn join("\t", $amino_acid_start, $amino_acid_end) . "\n";
			      foreach my $snps_entry (@{$snps_positions->{$query_id}}){
				    warn $snps_entry . "\n";
				    my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);
				    warn "($snp_position >= $query_start) and ($snp_position <= $query_end)" unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
				    next unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					
				    my @split_snp_variants = split(/\//, $snp_variants);
				    my ($snp_variant1, $snp_variant2) = ($split_snp_variants[0],$split_snp_variants[1]);
				    my $new_snp_variants = join("/", get_comp_base($snp_variant1), get_comp_base($snp_variant2));

				    my $new_snp_position = ($query_end - $snp_position) + $amino_acid_start;
# 				    die if ($new_snp_position < 0);
				    if(($new_snp_position >= $amino_acid_start) and ($new_snp_position <= $amino_acid_end)){
					  my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence) = get_aa_codon_snp($comp_query_subseq,$amino_acid_start,$new_snp_position,$new_snp_variants);
					  
					  my @snp_variants_list_entries = split(/\//, $snp_variants_list);
					  my $num_snp_variants = scalar(@snp_variants_list_entries);
					  if($num_snp_variants eq 3){
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"CDS_SNP",$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence"));
					  }else{

						$coding_region{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"CDS_SNP",$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence"));
					  }
					  
				    }elsif($new_snp_position < $amino_acid_start){
					  warn "5 PRIME UTR\n";
					  my $five_prime_UTR_start = $target_start;
					  my $five_prime_UTR_end = $target_start + ($amino_acid_start - 1);
      # 						my $five_prime_UTR_start = $target_start;
      # 						my $five_prime_UTR_end = $target_end;

					  my $antisense_seq = reverse_complement($query_seqs->{$query_id});
					  my $snp_base = get_subseq($antisense_seq, $new_snp_position, $new_snp_position);
					  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
						$five_prime_UTR{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"five_prime_UTR_SNP",$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$new_snp_variants]"));
					  }elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
						my $snp_variant3 = $snp_base;
						my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"five_prime_UTR_SNP",$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]"));
					  }
				    }
				    elsif($new_snp_position > $amino_acid_end){
					  warn "3 PRIME UTR\n";
					  my $three_prime_UTR_start = $target_start + ($amino_acid_end + 1);
					  my $three_prime_UTR_end = $target_end;
      # 						my $three_prime_UTR_start = $target_start;
      # 						my $three_prime_UTR_end = $target_end;
					  my $antisense_seq = reverse_complement($query_seqs->{$query_id});
					  my $snp_base = get_subseq($antisense_seq, $new_snp_position, $new_snp_position);
					  my ($snp_variant1, $snp_variant2) = split(/\//, $new_snp_variants);
					  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
						$three_prime_UTR{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"three_prime_UTR_SNP",($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$new_snp_variants]"));
					  }elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
						my $snp_variant3 = $snp_base;
						my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"three_prime_UTR_SNP",($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]"));
					  }
				    }
			      }
			}
			case "?" { 
			      warn "IN INDETERMINATE STRAND (?)\n";
# 			      warn $gmap_alignments->{$alignment_name} . "\n";
			      warn $query_id . "\n";
			      my @query_subseq = split('', get_subseq($query_seqs->{$query_id}, $amino_acid_start, $amino_acid_end));

			      my $query_subsequence = join('', @query_subseq);

			      my $aa_seq = translate_dna($query_subsequence);
# 			      warn "$query_id\n" . join('', @query_subseq) . "\n";
# 			      warn "$query_id\n$aa_seq\n";
			      warn join("\t", $amino_acid_start, $amino_acid_end) . "\n";
			      foreach my $snps_entry (@{$snps_positions->{$query_id}}){
				    warn $snps_entry . "\n";
				    my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);
				    warn "($snp_position >= $query_start) and ($snp_position <= $query_end)" unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
				    next unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					  
				    if(($snp_position >= $amino_acid_start) and ($snp_position <= $amino_acid_end)){
					  my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence) = get_aa_codon_snp($query_subsequence,$amino_acid_start,$snp_position,$snp_variants);
      
					  my @snp_variants_list_entries = split(/\//, $snp_variants_list);
					  my $num_snp_variants = scalar(@snp_variants_list_entries);
					  if($num_snp_variants eq 3){
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"CDS_SNP",$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence"));
					  }else{

						$coding_region{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"CDS_SNP",$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence"));
					  }
					  
				    }elsif($snp_position < $amino_acid_start){
					  warn "5 PRIME UTR\n";
					  my $five_prime_UTR_start = $target_start;
					  my $five_prime_UTR_end = $target_start + ($amino_acid_start - 1);
# 						my $five_prime_UTR_start = $target_start;
# 						my $five_prime_UTR_end = $target_end;

					  my $snp_base = get_subseq($query_seqs->{$query_id}, $snp_position, $snp_position);
					  my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
					  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
						$five_prime_UTR{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"five_prime_UTR_SNP",$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]"));
					  }elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
						my $snp_variant3 = $snp_base;
						my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"five_prime_UTR_SNP",$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]"));
					  }
				    }
				    elsif($snp_position > $amino_acid_end){
					  warn "3 PRIME UTR\n";
					  my $three_prime_UTR_start = $target_start + ($amino_acid_end + 1);
					  my $three_prime_UTR_end = $target_end;
# 						my $three_prime_UTR_start = $target_start;
# 						my $three_prime_UTR_end = $target_end;

					  my $snp_base = get_subseq($query_seqs->{$query_id}, $snp_position, $snp_position);
					  my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
					  if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
						$three_prime_UTR{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"three_prime_UTR_SNP",($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]"));
					  }elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
						my $snp_variant3 = $snp_base;
						my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
						$indeterminate{$alignment_name} = join("\t",$alignment_name,$coverage,$percent_identity,$query_id,"three_prime_UTR_SNP",($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]"));
					  }
				    }
	    # 			      
			      }
			}
			else { 
			      die "query_strand: $query_strand is not one of the following characters +, -, ?." 
			}
		  }
	    }
      }else{
	    warn "no snp information for $query_id";
      }
}


my $snp_filename = fileparse($gmap_infile, qr/\.[^.]*/);

my $snps_outfile = join('/', $output_dir, join("_", $snp_filename, "snps_outfile.txt"));
open(OUTFILE, ">$snps_outfile") or die "Couldn't open file $snps_outfile for writting, $!";                          
print OUTFILE join("\t", "alignment_name","coverage","percent_identity","query_id","snp_type","query_start","query_end","query_strand","target_id","target_start",
		  "target_end","target_strand","snp_description") . "\n";
foreach my $alignment_name (sort keys %coding_region){
      print OUTFILE $coding_region{$alignment_name} . "\n";
}
foreach my $alignment_name (sort keys %five_prime_UTR){
      print OUTFILE $five_prime_UTR{$alignment_name} . "\n";
}
foreach my $alignment_name (sort keys %three_prime_UTR){
      print OUTFILE $three_prime_UTR{$alignment_name} . "\n";
}
foreach my $alignment_name (sort keys %paralogs){
      print OUTFILE $paralogs{$alignment_name} . "\n";
}
close(OUTFILE) or die "Couldn't close file $snps_outfile";

my $indeterminate_outfile = join('/', $output_dir, join("_", $snp_filename, "indeterminate_outfile.txt"));
open(OUTFILE, ">$indeterminate_outfile") or die "Couldn't open file $indeterminate_outfile for writting, $!";                          
print OUTFILE join("\t", "alignment_name","coverage","percent_identity","query_id","snp_type","query_start","query_end","query_strand","target_id","target_start",
		  "target_end","target_strand","snp_description") . "\n";
foreach my $alignment_name (sort keys %indeterminate){
      print OUTFILE $indeterminate{$alignment_name} . "\n";
}
close(OUTFILE) or die "Couldn't close file $indeterminate_outfile";

sub parse_gmap_infile{
      my $gmap_infile = shift;
      die "Error lost gmap input file" unless defined $gmap_infile;

      my %gmap_alignments;
      my $i = 0;
      open(INFILE, "<$gmap_infile") or die "Couldn't open file $gmap_infile for reading, $!";
      while(<INFILE>){
	    chomp $_;
#  	    warn "$_\n";
	    if($i ne 0){
		  my @split_gmap_align_entry = split(/\t/, $_);
		  my ($alignment_name, $coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,$query_start,
			$query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,$target_strand,
			$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$query_align_block,
			$target_align_block,$align_identity_block,$intron_length_block) = @split_gmap_align_entry;

#  		  warn "$alignment_name\n";

		  my $gmap_align_entry = join("\t", $coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,
			$query_start,$query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,
			$target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,
			$query_align_block,$target_align_block,$align_identity_block,$intron_length_block);

		  $gmap_alignments{$alignment_name} = $gmap_align_entry;
	    }
	    $i++;
      }
      close(INFILE) or die "Couldn't close file $gmap_infile";
      
      return \%gmap_alignments;
}

sub parse_snps_infile{
      my $snps_infile = shift;
      die "Error lost gmap input file" unless defined $snps_infile;

      my %snp_positions;
      my $i = 0;
      open(INFILE, "<$snps_infile") or die "Couldn't open file $snps_infile for reading, $!";
      while(<INFILE>){
	    chomp $_;
#  	    warn "$_\n";
	    if($i ne 0){
		  my @split_entry = split(/\t/, $_);
		  my ($query_id, $snp_position, $snp_variants) = ($split_entry[0], $split_entry[1], $split_entry[2]); 
# 		  warn "$query_id\n";
		  my $snp_entry = join("\t", $snp_position, $snp_variants);
		  push(@{$snp_positions{$query_id}}, $snp_entry);
	    }
	    $i++;
      }
      close(INFILE) or die "Couldn't close file $snps_infile";
      
      return \%snp_positions;
}


sub parse_fasta_infile{

      my $fasta_infile = shift;
      die "Error lost fasta input file" unless defined $fasta_infile;

      my %fasta_seqs;
      my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
      while(my $seq_entry = $seqio->next_seq) {

	    my $seq_id = $seq_entry->id;
	    my $sequence = $seq_entry->seq;

# 	    warn $seq_id . "\n";
# 	    warn $sequence . "\n";
	    $fasta_seqs{$seq_id} = $sequence;
      }

      # Clean out the sequence I/O object.
      $seqio = ();

      return \%fasta_seqs;
}

# my $seq = get_subseq("AGCTTGCGTT", 3, 8);
# warn $seq . "\n";
sub get_subseq{

      my $sequence = shift;
      die "Error lost sequence" unless defined $sequence;

      my $seq_start = shift;
      die "Error lost start of sequence" unless defined $seq_start;

      my $seq_end = shift;
      die "Error lost end of sequence" unless defined $seq_end;

      $seq_start = $seq_start - 1;
      $seq_end = $seq_end;

      my $length = ($seq_end - $seq_start);

      my $trimmed_seq = substr($sequence, $seq_start, $length);
      
      return uc($trimmed_seq);
}

# my $protein_sequence = translate_dna("GTGGCTCCAGGGAATGCACTTGCAATTGTTAATGCAACTGCAGGGTCTACAACTGTGGGAGCAGAAAGGGGTAGTGCTTCAGTTTCCATCTCTGGGCTT");
# warn "VAPGNALAIVNATAGSTTVGAERGSASVSISGL\n";
# warn $protein_sequence . "\n";

sub translate_dna{

      my $dna_sequence = shift;
      die "Error lost input dna sequence" unless defined $dna_sequence;

      my $aa_sequence = "";
      for(my $i = 0; $i < (length($dna_sequence) - 2); $i += 3) {
	    my $codon = substr($dna_sequence, $i, 3);
	    $aa_sequence .= convert_codon($codon);
      }

      return $aa_sequence;
}

sub get_aa_codon_snp{

      my $dna_sequence = shift;
      die "Error lost input dna sequence" unless defined $dna_sequence;
      
      my $amino_acid_start = shift;
      die "Error lost input amino acid start position" unless defined $amino_acid_start;
      
      my $snp_position = shift;
      die "Error lost input snp position" unless defined $snp_position;

      my $snp_variants = shift;
      die "Error lost input snp variants" unless defined $snp_variants;

      my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);

      my $adj_snp_position = $snp_position - $amino_acid_start;

      my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $aa_sequence);

      my $subsitution_position = 0;
      my $is_snp_found = "false";

      my $amino_acid_position = 0;
      for(my $i = 0; $i < (length($dna_sequence) - 2); $i += 3) {
	    if(($i > $adj_snp_position) and ($is_snp_found eq "false")){
		  my $snp_codon_position = ($snp_position - (($amino_acid_start + $i) - 3));

		  warn join("\t", $snp_position,(($amino_acid_start + $i) - 3), $snp_codon_position) . "\n";

		  my $snp_codon = substr($dna_sequence, ($i - 3), 3);
		  my @snp_codon_bases = split('', $snp_codon);
		  
		  warn "snp_variant1 $snp_variant1 eq $snp_codon_bases[$snp_codon_position]\n";
		  warn "snp_variant2 $snp_variant2 eq $snp_codon_bases[$snp_codon_position]\n";
		  my ($snp_codon_variant1, $snp_codon_variant2, $snp_codon_variant3) = "";
		  my ($amino_acid1, $amino_acid2, $amino_acid3) = "";

		  if($snp_variant1 eq $snp_codon_bases[$snp_codon_position]){
 			warn "$snp_variant1 eq $snp_codon_bases[$snp_codon_position]\n";
			$snp_codon_variant1 = $snp_codon;
			$snp_codon_bases[$snp_codon_position] = $snp_variant2;
			$snp_codon_variant2 = join('', @snp_codon_bases);
			($amino_acid1, $amino_acid2) = (convert_codon($snp_codon_variant1), convert_codon($snp_codon_variant2));
			if($amino_acid1 eq $amino_acid2){
			      $substitution_type = "synonymous";
			}else{
			      $substitution_type = "nonsynonymous";
			}
			$codon_list = join("/", $snp_codon_variant1, $snp_codon_variant2);
			$amino_acid_list = join("/", $amino_acid1, $amino_acid2);
			my @split_amino_acid_list = split(/\//, $amino_acid_list);
			my @converted_aa_bases = ();
			foreach my $aa_base (@split_amino_acid_list){
			      push(@converted_aa_bases, convert_aa_base($aa_base));
			}
			$amino_acid_list = join("/", @converted_aa_bases);
			
		  }elsif($snp_variant2 eq $snp_codon_bases[$snp_codon_position]){
 			warn "$snp_variant2 eq $snp_codon_bases[$snp_codon_position]\n";
			$snp_codon_variant1 = $snp_codon;
			$snp_codon_bases[$snp_codon_position] = $snp_variant1;
			$snp_codon_variant2 = join('', @snp_codon_bases);
			($amino_acid1, $amino_acid2) = (convert_codon($snp_codon_variant1), convert_codon($snp_codon_variant2));
			if($amino_acid1 eq $amino_acid2){
			      $substitution_type = "synonymous";
			}else{
			      $substitution_type = "nonsynonymous";
			}
			$codon_list = join("/", $snp_codon_variant1, $snp_codon_variant2);
			$amino_acid_list = join("/", $amino_acid1, $amino_acid2);
			my @split_amino_acid_list = split(/\//, $amino_acid_list);
			my @converted_aa_bases = ();
			foreach my $aa_base (@split_amino_acid_list){
			      push(@converted_aa_bases, convert_aa_base($aa_base));
			}
			$amino_acid_list = join("/", @converted_aa_bases);

		  }else{ # complex/indeterminate snp
			
			my $snp_variant3 = $snp_codon_bases[$snp_codon_position];
			$snp_codon_variant3 = $snp_codon;
			$snp_codon_bases[$snp_codon_position] = $snp_variant1;
			$snp_codon_variant1 = join('', @snp_codon_bases);
			$snp_codon_bases[$snp_codon_position] = $snp_variant2;
			$snp_codon_variant2 = join('', @snp_codon_bases);
			($amino_acid1, $amino_acid2, $amino_acid3) = (convert_codon($snp_codon_variant1), convert_codon($snp_codon_variant2), convert_codon($snp_codon_variant3));

			if(($amino_acid1 eq $amino_acid3) or ($amino_acid2 eq $amino_acid3)){
			      $substitution_type = "synonymous";
			}else{
			      $substitution_type = "nonsynonymous";
			}
			$codon_list = join("/", $snp_codon_variant1, $snp_codon_variant2, $snp_codon_variant3);
			$amino_acid_list = join("/", $amino_acid1, $amino_acid2, $amino_acid3);
			my @split_amino_acid_list = split(/\//, $amino_acid_list);
			my @converted_aa_bases = ();
			foreach my $aa_base (@split_amino_acid_list){
			      push(@converted_aa_bases, convert_aa_base($aa_base));
			}
			$amino_acid_list = join("/", @converted_aa_bases);
			$snp_variants_list = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
		  }
		  
		  $is_snp_found = "true";
		  
		  # The substitution position is shifted to the left in this if statement.
		  $subsitution_position = ($amino_acid_position - 1);
	    }
	    my $codon = substr($dna_sequence, $i, 3);
	    $aa_sequence .= convert_codon($codon);
	    $amino_acid_position++;
      }
      my @split_aa_sequence = split('', $aa_sequence);
      my $amino_acid_base = $split_aa_sequence[$subsitution_position];
      $split_aa_sequence[$subsitution_position] = "[$amino_acid_base]";
      my $amino_acid_sequence = join('', @split_aa_sequence);

      $snp_variants_list = $snp_variants unless(defined($snp_variants_list));
      return ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence);
}


sub convert_codon{
      my $codon = shift;
      die "Error lost codon" unless defined $codon;
      
      my $amino_acid = "";
      switch ($codon) {
	    case /GCT|GCC|GCA|GCG/		{ $amino_acid = "A" }
	    case /TGT|TGC/			{ $amino_acid = "C" }
	    case /GAT|GAC/			{ $amino_acid = "D" }
	    case /GAA|GAG/			{ $amino_acid = "E" }
	    case /TTT|TTC/			{ $amino_acid = "F" }
	    case /GGT|GGC|GGA|GGG/		{ $amino_acid = "G" }
	    case /CAT|CAC/			{ $amino_acid = "H" }
	    case /ATT|ATC|ATA/			{ $amino_acid = "I" }
	    case /AAA|AAG/			{ $amino_acid = "K" }
	    case /TTA|TTG|CTT|CTC|CTA|CTG/	{ $amino_acid = "L" }
	    case /ATG/				{ $amino_acid = "M" }
	    case /AAT|AAC/			{ $amino_acid = "N" }
	    case /CCT|CCC|CCA|CCG/		{ $amino_acid = "P" }
	    case /CAA|CAG/			{ $amino_acid = "Q" }
	    case /CGT|CGC|CGA|CGG|AGA|AGG/	{ $amino_acid = "R" }
	    case /TCT|TCC|TCA|TCG|AGT|AGC/	{ $amino_acid = "S" }
	    case /ACT|ACC|ACA|ACG/		{ $amino_acid = "T" }
	    case /GTT|GTC|GTA|GTG/		{ $amino_acid = "V" }
	    case /TGG/				{ $amino_acid = "W" }
	    case /TAT|TAC/			{ $amino_acid = "Y" }
	    case /TAA|TGA|TAG/			{ $amino_acid = "*" }
	    else 				{ $amino_acid = "X" }
      }
      return $amino_acid;
}

sub convert_aa_base{

      my $aa_base = shift;
      die "Error lost amino acid base" unless defined $aa_base;
      
      my $aa_abbrev = "";
      switch ($aa_base) {
	    case /A/	{ $aa_abbrev = "Ala" }
	    case /R/	{ $aa_abbrev = "Arg" }
	    case /N/	{ $aa_abbrev = "Asn" }
	    case /D/	{ $aa_abbrev = "Asp" }
	    case /C/	{ $aa_abbrev = "Cys" }
	    case /E/	{ $aa_abbrev = "Glu" }
	    case /Q/	{ $aa_abbrev = "Gln" }
	    case /G/	{ $aa_abbrev = "Gly" }
	    case /H/	{ $aa_abbrev = "His" }
	    case /I/	{ $aa_abbrev = "Ile" }
	    case /L/	{ $aa_abbrev = "Leu" }
	    case /K/	{ $aa_abbrev = "Lys" }
	    case /M/	{ $aa_abbrev = "Met" }
	    case /F/	{ $aa_abbrev = "Phe" }
	    case /P/	{ $aa_abbrev = "Pro" }
	    case /S/	{ $aa_abbrev = "Ser" }
	    case /T/	{ $aa_abbrev = "Thr" }
	    case /W/	{ $aa_abbrev = "Trp" }
	    case /Y/	{ $aa_abbrev = "Tyr" }
	    case /V/	{ $aa_abbrev = "Val" }
	    case /X/	{ $aa_abbrev = "Unknown" }
	    case "*"	{ $aa_abbrev = "Stop" }
	    else { die "aa_base: $aa_base is not one of the following characters ACDEFGHIKLMNPQRSTVWXY*"  }	    
      }
      return $aa_abbrev;
}

sub reverse_complement{

      my $dna_sequence = shift;
      die "Error lost input dna sequence" unless defined $dna_sequence;

      my @comp_seq = ();
      my @dna_seq = split('', $dna_sequence);
      for(my $i = 0; $i < length($dna_sequence); $i++) {
	    my $nuc_base = $dna_seq[$i];
	    my $comp_base = get_comp_base($nuc_base);
	    push(@comp_seq, $comp_base);
      }


      return reverse(@comp_seq);
}

sub get_comp_base{
      my $nuc_base = shift;
      die "Error lost input nucleotide base" unless defined $nuc_base;

      my $comp_base = "";
      switch ($nuc_base) {
	    case /A/		{ $comp_base = "T" }
	    case /C/		{ $comp_base = "G" }
	    case /G/		{ $comp_base = "C" }
	    case /T/		{ $comp_base = "A" }
	    case m/[RYKMSWBDHVN]/i {$comp_base = $nuc_base}
	    else { die "nuc_base: $nuc_base is not one of the following characters ACGTRYKMSWBDHVN"  }
      }
      return $comp_base;
}