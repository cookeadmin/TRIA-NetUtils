#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use IPC::Open2;
use Bio::SeqIO;
use File::Basename;
use List::Compare;
use List::Util qw(min max);
use Switch;

# perl gmap_snps_analysis.pl -s /home/cookeadmin/workspace/cathy/AllSNPsVariantConsensusPositions.txt -g /home/cookeadmin/workspace/cathy/AllSNPs_output_GMAP.parsed -q /home/cookeadmin/workspace/cathy/AllSNPsContigs.fa -o ~/workspace/cathy
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


my ($query_seqs, $query_lengths) = parse_fasta_infile($query_infile);

my $gmap_alignments = parse_gmap_infile($gmap_infile);

my $snps_positions = parse_snps_infile($snps_infile);

my (%gmap_snps_output, @snps_positions_list, @gmap_snps_list, @snps_gmap_alignment, @snps_no_gmap_alignment, @gmap_no_snps) = ();

# The entire SNP positions list to find the difference between the GMAP SNP positions list and the SNP positions list.
@snps_positions_list = keys %{$snps_positions};

foreach my $alignment_name (sort keys %{$gmap_alignments}){
	warn $alignment_name . "\n";
	my @split_gmap_align_entry = split(/\t/, $gmap_alignments->{$alignment_name});
	my ($coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,$query_start,
		$query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,$target_strand,
		$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$query_align_block,
		$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,
		$query_gap_position_list,$mismatch_position_list) = @split_gmap_align_entry;

	warn join("\t", $alignment_name,$coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,$query_start,
		$query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,$target_strand,
		$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$query_align_block,
		$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,
		$query_gap_position_list,$mismatch_position_list) . "\n";
	
	# If we have any gaps within the query sequence we need to figure out if we need adjust the amino acid start or amino acid end position to fix the possibility of a frame shift 
	# causing the amino acid sequence to be incorrect. 
	# GMAP only evaluates where the amino acid translation starts in the target sequence, which may introduce a frame shift.
	my $query_gap_amino_acid_range_counter = 0;
	if($query_gap_position_list ne "N/A"){
		my @query_gap_positions = split(/,/, $query_gap_position_list);
		foreach my $query_gap_position (@query_gap_positions){
			$query_gap_amino_acid_range_counter++ if(($query_gap_position >= $amino_acid_start) and ($query_gap_position <= $amino_acid_end));
		}
	}

	# If we have any gaps within the target sequence we need to figure out if we need adjust the amino acid start or amino acid end position to fix the possibility of a frame shift 
	# causing the amino acid sequence to be incorrect. 
	# GMAP only evaluates where the amino acid translation starts in the target sequence, which may introduce a frame shift.
	my $target_gap_amino_acid_range_counter = 0;
	if($target_gap_position_list ne "N/A"){
		my @target_gap_positions = split(/,/, $target_gap_position_list);
		foreach my $target_gap_position (@target_gap_positions){
			$target_gap_amino_acid_range_counter++ if(($target_gap_position >= $amino_acid_start) and ($target_gap_position <= $amino_acid_end));
		}
	}
# 	die $target_gap_amino_acid_range_counter if(($query_strand eq "-") and ($target_gap_amino_acid_range_counter >= 1) and ($query_id eq "JackDenovo_rep_c1237"));
# 		die join("; ", $query_id, $target_gap_amino_acid_range_counter, ($target_gap_amino_acid_range_counter % 3)) if(($query_strand eq "-") and (($target_gap_amino_acid_range_counter % 3) ne 0));

	
	
	# Get the number of alignment paths and the current path number.
	my ($path_num, $num_align_paths) = "";
	if($alignment_name =~ m/path(\d+)of(\d+)/){
		($path_num, $num_align_paths) = ($1, $2);
		warn join("\t", $path_num, $num_align_paths) . "\n";
	}

	# Decide whether or not the SNP is a possible paralog based on more than one path. 
	my $snp_type = "SNP";
	$snp_type = "SNP_Paralog" if ($num_align_paths > 1);

	# If the number of alignment paths are greater than one and the position of the SNP information for the query id exists.
	if(($num_align_paths >= 1) and defined(@{$snps_positions->{$query_id}})){
		# The entire GMAP SNP positions list to find the difference between the SNP positions list and the GMAP SNP positions list.
		push(@gmap_snps_list, $query_id);
		switch ($query_strand) {
			case "+" {
				warn "IN SENSE STRAND (+)\n";
		# 		warn $gmap_alignments->{$alignment_name} . "\n";
				warn $query_id . "\n";
				my @query_subseq = split('', get_subseq($query_seqs->{$query_id}, $amino_acid_start, $amino_acid_end));

				my $query_subsequence = join('', @query_subseq);

				my $aa_seq = translate_dna($query_subsequence);
		# 		warn "$query_id\n" . join('', @query_subseq) . "\n";
		# 		warn "$query_id\n$aa_seq\n";
		# 		warn join("\t", $amino_acid_start, $amino_acid_end) . "\n";
				foreach my $snps_entry (sort @{$snps_positions->{$query_id}}){
					warn $snps_entry . "\n";
					my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);

					warn "($snp_position >= $query_start) and ($snp_position <= $query_end)" unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					next unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					
					# If the amino acid translation information exists process as 5' or 3' untranslated region or coding region based on the position of the SNP.
					if(($amino_acid_start ne "N/A") and ($amino_acid_end ne "N/A") and ($amino_acid_length ne "N/A")){
						if(($snp_position >= $amino_acid_start) and ($snp_position <= $amino_acid_end)){
							my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence) = get_aa_codon_snp($query_subsequence,$query_strand,$amino_acid_start,$snp_position,$snp_variants);
							
							my @snp_variants_list_entries = split(/\//, $snp_variants_list);
							my $num_snp_variants = scalar(@snp_variants_list_entries);
							if($num_snp_variants eq 3){
								push(@{$gmap_snps_output{$alignment_name}{"CDS"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "CDS","Complex", $snp_type),$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence")));
							}else{
								push(@{$gmap_snps_output{$alignment_name}{"CDS"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "CDS", $snp_type),$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence")));
							}
						}elsif($snp_position < $amino_acid_start){
							warn "5 PRIME UTR\n";
							my $five_prime_UTR_start = $target_start;
							my $five_prime_UTR_end = $target_start + ($amino_acid_start - 1);
	# 						my $five_prime_UTR_start = $target_start;
	# 						my $five_prime_UTR_end = $target_end;

							my $query_sequence = $query_seqs->{$query_id};
							my @query_seq = split('', $query_sequence);
							my $snp_base = $query_seq[$snp_position - 1];

							my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
							if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
								push(@{$gmap_snps_output{$alignment_name}{"five_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "five_prime_UTR", $snp_type),$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
							}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
								my $snp_variant3 = $snp_base;
								my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
								push(@{$gmap_snps_output{$alignment_name}{"five_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "five_prime_UTR", "Complex", $snp_type),$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
							}
						}
						elsif($snp_position > $amino_acid_end){
							warn "3 PRIME UTR\n";
							my $three_prime_UTR_start = $target_start + ($amino_acid_end + 1);
							my $three_prime_UTR_end = $target_end;
	# 						my $three_prime_UTR_start = $target_start;
	# 						my $three_prime_UTR_end = $target_end;

							my $query_sequence = $query_seqs->{$query_id};
							my @query_seq = split('', $query_sequence);
							my $snp_base = $query_seq[$snp_position - 1];

							my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
							if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
								push(@{$gmap_snps_output{$alignment_name}{"three_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "three_prime_UTR", $snp_type),($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
							}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
								my $snp_variant3 = $snp_base;
								my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
								push(@{$gmap_snps_output{$alignment_name}{"three_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "three_prime_UTR", "Complex", $snp_type),($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
							}
						}
					}else{ 
						warn "NONCODING\n";

						my $query_sequence = $query_seqs->{$query_id};
						my @query_seq = split('', $query_sequence);
						my $snp_base = $query_seq[$snp_position - 1];
						my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
						if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
							push(@{$gmap_snps_output{$alignment_name}{"Noncoding"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "Noncoding", $snp_type),$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
						}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
							my $snp_variant3 = $snp_base;
							my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
							push(@{$gmap_snps_output{$alignment_name}{"Noncoding"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "Noncoding", "Complex", $snp_type),$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
						}
					}
				}
			}
			case "-" { 
				warn "IN ANTISENSE STRAND (-)\n";
		# 		warn $gmap_alignments->{$alignment_name} . "\n";
				warn $query_id . "\n";
# 				if((($target_gap_amino_acid_range_counter % 3) ne 0) and (($query_gap_amino_acid_range_counter % 3) eq 0)){
# 					$amino_acid_start = ($amino_acid_start + ($target_gap_amino_acid_range_counter % 3));
# 				}elsif((($target_gap_amino_acid_range_counter % 3) eq 0) and (($query_gap_amino_acid_range_counter % 3) ne 0)){
# 					$amino_acid_start = ($amino_acid_start - ($query_gap_amino_acid_range_counter % 3));
# 				}els
				if((($target_gap_amino_acid_range_counter % 3) > ($query_gap_amino_acid_range_counter % 3))){
					$amino_acid_start = ($amino_acid_start + (($target_gap_amino_acid_range_counter % 3) - ($query_gap_amino_acid_range_counter % 3)));
				}elsif((($target_gap_amino_acid_range_counter % 3) < ($query_gap_amino_acid_range_counter % 3))){
					$amino_acid_start = ($amino_acid_start - (($query_gap_amino_acid_range_counter % 3) - ($target_gap_amino_acid_range_counter % 3))) if($amino_acid_start > 2);
					$amino_acid_start = ($amino_acid_start + (($query_gap_amino_acid_range_counter % 3) - ($target_gap_amino_acid_range_counter % 3))) if($amino_acid_start <= 2);
				}
# 				warn join("; ", $query_id, $target_gap_amino_acid_range_counter, ($target_gap_amino_acid_range_counter % 3), $query_gap_amino_acid_range_counter, ($query_gap_amino_acid_range_counter % 3)) if($query_strand eq "-");
# die if $query_id eq "LodgeDenovo_rep_c103";
				my @query_subseq = split('', get_subseq($query_seqs->{$query_id}, $amino_acid_start, $amino_acid_end));

				my $query_subsequence = join('', @query_subseq);

				my $aa_seq = translate_dna($query_subsequence);
		# 		warn "$query_id\n" . join('', @query_subseq) . "\n";
		# 		warn "$query_id\n$aa_seq\n";
		# 		warn join("\t", $amino_acid_start, $amino_acid_end) . "\n";
				foreach my $snps_entry (sort @{$snps_positions->{$query_id}}){
					warn $snps_entry . "\n";
					my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);

					warn "($snp_position >= $query_start) and ($snp_position <= $query_end)" unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					next unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					# If the amino acid translation information exists process as 5' or 3' untranslated region or coding region based on the position of the SNP.
					if(($amino_acid_start ne "N/A") and ($amino_acid_end ne "N/A") and ($amino_acid_length ne "N/A")){
						if(($snp_position >= $amino_acid_start) and ($snp_position <= $amino_acid_end)){
							my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence) = get_aa_codon_snp($query_subsequence,$query_strand,$amino_acid_start,$snp_position,$snp_variants);
							my @snp_variants_list_entries = split(/\//, $snp_variants_list);
							my $num_snp_variants = scalar(@snp_variants_list_entries);
							if($num_snp_variants eq 3){
								push(@{$gmap_snps_output{$alignment_name}{"CDS"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "CDS", "Complex", $snp_type),$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence")));
							}else{
								push(@{$gmap_snps_output{$alignment_name}{"CDS"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "CDS", $snp_type),$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence")));
							}
							
						}
						elsif($snp_position > $amino_acid_end){
							warn "5 PRIME UTR\n";

							my $five_prime_UTR_start = $target_end;
							my $five_prime_UTR_end = $target_end - ($amino_acid_end + 1);
			# 				my $five_prime_UTR_start = $target_start;
			# 				my $five_prime_UTR_end = $target_end;
							my $query_sequence = $query_seqs->{$query_id};
							$query_sequence = join('', complement($query_sequence)) if($query_strand eq "-");
							my @query_seq = split('', $query_sequence);
							my $snp_base = $query_seq[$snp_position - 1];

							my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
							($snp_variant1, $snp_variant2) = (get_comp_base($snp_variant1), get_comp_base($snp_variant2)) if($query_strand eq "-");
							$snp_variants = join("/", $snp_variant1, $snp_variant2) if($query_strand eq "-");
							if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
								push(@{$gmap_snps_output{$alignment_name}{"five_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "five_prime_UTR", $snp_type),($amino_acid_end + 1),$query_end,$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
							}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
								my $snp_variant3 = $snp_base;
								my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
								push(@{$gmap_snps_output{$alignment_name}{"five_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "five_prime_UTR", "Complex", $snp_type),($amino_acid_end + 1),$query_end,$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
							}
						}
						elsif($snp_position < $amino_acid_start){
							warn "3 PRIME UTR\n";

							my $three_prime_UTR_start = $target_start;
							my $three_prime_UTR_end = $target_start + ($amino_acid_start - 1);
			# 				my $three_prime_UTR_start = $target_start;
			# 				my $three_prime_UTR_end = $target_end;
							my $query_sequence = $query_seqs->{$query_id};
							$query_sequence = join('', complement($query_sequence)) if($query_strand eq "-");
							my @query_seq = split('', $query_sequence);
							my $snp_base = $query_seq[$snp_position - 1];

							my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
							($snp_variant1, $snp_variant2) = (get_comp_base($snp_variant1), get_comp_base($snp_variant2)) if($query_strand eq "-");
							$snp_variants = join("/", $snp_variant1, $snp_variant2) if($query_strand eq "-");
							if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
								push(@{$gmap_snps_output{$alignment_name}{"three_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "three_prime_UTR", $snp_type),$query_start,($amino_acid_start - 1),$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
							}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
								my $snp_variant3 = $snp_base;
								my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
								push(@{$gmap_snps_output{$alignment_name}{"three_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "three_prime_UTR", "Complex", $snp_type),$query_start,($amino_acid_start - 1),$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
							}
						}
					}else{ 
						warn "NONCODING\n";

						my $query_sequence = $query_seqs->{$query_id};
						my @query_seq = split('', $query_sequence);
						my $snp_base = $query_seq[$snp_position - 1];
						my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
						if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
							push(@{$gmap_snps_output{$alignment_name}{"Noncoding"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "Noncoding", $snp_type),$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
						}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
							my $snp_variant3 = $snp_base;
							my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
							push(@{$gmap_snps_output{$alignment_name}{"Noncoding"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "Noncoding", "Complex", $snp_type),$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
						}
					}
				}
			}
			case "?" { 
				warn "IN INDETERMINATE STRAND (?)\n";
		# 		warn $gmap_alignments->{$alignment_name} . "\n";
				warn $query_id . "\n";
				my @query_subseq = split('', get_subseq($query_seqs->{$query_id}, $amino_acid_start, $amino_acid_end));

				my $query_subsequence = join('', @query_subseq);

				my $aa_seq = translate_dna($query_subsequence);
		# 		warn "$query_id\n" . join('', @query_subseq) . "\n";
		# 		warn "$query_id\n$aa_seq\n";
				warn join("\t", $amino_acid_start, $amino_acid_end) . "\n";
				foreach my $snps_entry (sort @{$snps_positions->{$query_id}}){
					warn $snps_entry . "\n";
					my ($snp_position, $snp_variants) = split(/\t/, $snps_entry);
					warn "($snp_position >= $query_start) and ($snp_position <= $query_end)" unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					next unless (($snp_position >= $query_start) and ($snp_position <= $query_end));
					# If the amino acid translation information exists process as 5' or 3' untranslated region or coding region based on the position of the SNP.
					if(($amino_acid_start ne "N/A") and ($amino_acid_end ne "N/A") and ($amino_acid_length ne "N/A")){
						if(($snp_position >= $amino_acid_start) and ($snp_position <= $amino_acid_end)){
							my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $amino_acid_sequence) = get_aa_codon_snp($query_subsequence,$query_strand,$amino_acid_start,$snp_position,$snp_variants);

							my @snp_variants_list_entries = split(/\//, $snp_variants_list);
							my $num_snp_variants = scalar(@snp_variants_list_entries);
							if($num_snp_variants eq 3){
								push(@{$gmap_snps_output{$alignment_name}{"CDS"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "CDS","Complex", $snp_type),$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence")));
							}else{
								push(@{$gmap_snps_output{$alignment_name}{"CDS"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "CDS", $snp_type),$amino_acid_start,$amino_acid_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants_list]","codons=[$codon_list]","substitution=[$amino_acid_list]","type=$substitution_type","sequence=$amino_acid_sequence")));
							}
						}elsif($snp_position < $amino_acid_start){
							warn "5 PRIME UTR\n";
							my $five_prime_UTR_start = $target_start;
							my $five_prime_UTR_end = $target_start + ($amino_acid_start - 1);
		# 					my $five_prime_UTR_start = $target_start;
		# 					my $five_prime_UTR_end = $target_end;

							my $query_sequence = $query_seqs->{$query_id};
							my @query_seq = split('', $query_sequence);
							my $snp_base = $query_seq[$snp_position - 1];
							my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
							if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
								push(@{$gmap_snps_output{$alignment_name}{"five_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "five_prime_UTR", $snp_type),$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
							}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
								my $snp_variant3 = $snp_base;
								my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
								push(@{$gmap_snps_output{$alignment_name}{"five_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "five_prime_UTR", "Complex", $snp_type),$query_start,($amino_acid_start - 1),$query_strand,$target_id,$five_prime_UTR_start,$five_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
							}
						}
						elsif($snp_position > $amino_acid_end){
							warn "3 PRIME UTR\n";
							my $three_prime_UTR_start = $target_start + ($amino_acid_end + 1);
							my $three_prime_UTR_end = $target_end;
		# 					my $three_prime_UTR_start = $target_start;
		# 					my $three_prime_UTR_end = $target_end;

							my $query_sequence = $query_seqs->{$query_id};
							my @query_seq = split('', $query_sequence);
							my $snp_base = $query_seq[$snp_position - 1];
							my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
							if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
								push(@{$gmap_snps_output{$alignment_name}{"three_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "three_prime_UTR", $snp_type),($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
							}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
								my $snp_variant3 = $snp_base;
								my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
								push(@{$gmap_snps_output{$alignment_name}{"three_prime_UTR"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "three_prime_UTR", "Complex", $snp_type),($amino_acid_end + 1),$query_end,$query_strand,$target_id,$three_prime_UTR_start,$three_prime_UTR_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
							}
						}
					}else{ 
						warn "NONCODING\n";

						my $query_sequence = $query_seqs->{$query_id};
						my @query_seq = split('', $query_sequence);
						my $snp_base = $query_seq[$snp_position - 1];
						my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
						if(($snp_base eq $snp_variant1) or ($snp_base eq $snp_variant2)){
							push(@{$gmap_snps_output{$alignment_name}{"Noncoding"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "Noncoding", $snp_type),$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$snp_variants]")));
						}elsif(($snp_base ne $snp_variant1) and ($snp_base ne $snp_variant2)){
							my $snp_variant3 = $snp_base;
							my $complex_snp_variants = join("/", $snp_variant1, $snp_variant2, $snp_variant3);
							push(@{$gmap_snps_output{$alignment_name}{"Noncoding"}}, join("\t",$alignment_name,$coverage,$percent_identity,$query_id,join("_", "Noncoding", "Complex", $snp_type),$query_start,$query_end,$query_strand,$target_id,$target_start,$target_end,$target_strand,join(";", "position=$snp_position","variants=[$complex_snp_variants]")));
						}
					}
		# 			      
				}
			}else{ 
				die "query_strand: $query_strand is not one of the following characters +, -, ?." 
			}
		}
	}else{
		warn "no snp information for $query_id";
		push(@gmap_no_snps, $query_id);
	}
}

my @gene_structures = ("five_prime_UTR", "CDS", "three_prime_UTR", "Noncoding");
my ($snp_counter, $paralog_snp_counter, $complex_snp_counter, $five_prime_UTR_snp_counter, $CDS_snp_counter, $three_prime_UTR_snp_counter, $noncoding_snp_counter) = 0;
my $snp_filename = fileparse($gmap_infile, qr/\.[^.]*/);

my $snps_outfile = join('/', $output_dir, join("_", $snp_filename, "snps_outfile.txt"));
warn "Writting the find GMAP snps outfile to $snps_outfile.....\n";
open(OUTFILE, ">$snps_outfile") or die "Couldn't open file $snps_outfile for writting, $!";                          
print OUTFILE join("\t", "alignment_name","coverage","percent_identity","query_id","snp_type","query_start","query_end","query_strand","target_id","target_start",
		"target_end","target_strand","snp_description") . "\n";
foreach my $alignment_name (sort keys %gmap_snps_output){
	foreach my $gene_structure (@gene_structures){
		foreach my $gmap_snps_output_entry (@{$gmap_snps_output{$alignment_name}{$gene_structure}}){
			if(defined($gmap_snps_output_entry)){
				print OUTFILE $gmap_snps_output_entry . "\n";
				$paralog_snp_counter++ if($gmap_snps_output_entry =~ m/SNP_Paralog/);
				$complex_snp_counter++ if($gmap_snps_output_entry =~ m/Complex_SNP/);
				$five_prime_UTR_snp_counter++ if($gene_structure eq "five_prime_UTR");
				$CDS_snp_counter++ if($gene_structure eq "CDS");
				$three_prime_UTR_snp_counter++ if($gene_structure eq "three_prime_UTR");
				$noncoding_snp_counter++ if($gene_structure eq "Noncoding");
				$snp_counter++;
			}
		}
	}
}
close(OUTFILE) or die "Couldn't close file $snps_outfile";


my $snps_stats_outfile = join('/', $output_dir, join("_", $snp_filename, "snps_stats_outfile.txt"));
warn "Writting the find GMAP SNPs Analysis Stats outfile to $snps_stats_outfile.....\n";
open(OUTFILE, ">$snps_stats_outfile") or die "Couldn't open file $snps_stats_outfile for writting, $!";

# Find all SNPs that do not have a GMAP alignment associated with them.
my %seen = ();
@snps_positions_list = grep { ! $seen{$_}++ } @snps_positions_list;

%seen = ();
@gmap_snps_list = grep { ! $seen{$_}++ } @gmap_snps_list;

my %list_comparison = ();
$list_comparison{$_}++ for(@snps_positions_list, @gmap_snps_list);

foreach my $entry (keys %list_comparison){
	push(@snps_no_gmap_alignment, $entry) if ($list_comparison{$entry} == 1);
	push(@snps_gmap_alignment, $entry) if ($list_comparison{$entry} == 2);
}

print OUTFILE join(" ", "Total number of processed SNPs:", $snp_counter) . "\n";
print OUTFILE join(" ", "5' UTR SNPs:", $five_prime_UTR_snp_counter, "out of", $snp_counter, "(" . sprintf("%.2f", ($five_prime_UTR_snp_counter/$snp_counter) * 100) . "%)") . "\n";
print OUTFILE join(" ", "Coding DNA Sequence SNPs:", $CDS_snp_counter, "out of", $snp_counter, "(" . sprintf("%.2f", ($CDS_snp_counter/$snp_counter) * 100) . "%)") . "\n";
print OUTFILE join(" ", "3' UTR SNPs:", $three_prime_UTR_snp_counter, "out of", $snp_counter, "(" . sprintf("%.2f", ($three_prime_UTR_snp_counter/$snp_counter) * 100) . "%)") . "\n";
print OUTFILE join(" ", "Noncoding SNPs:", $noncoding_snp_counter, "out of", $snp_counter, "(" . sprintf("%.2f", ($noncoding_snp_counter/$snp_counter) * 100) . "%)") . "\n";
print OUTFILE join(" ", "Possible Paralogous SNPs:", $paralog_snp_counter, "out of", $snp_counter, "(" . sprintf("%.2f", ($paralog_snp_counter/$snp_counter) * 100) . "%)") . "\n";
print OUTFILE join(" ", "Complex SNPs:", $complex_snp_counter, "out of", $snp_counter, "(" . sprintf("%.2f", ($complex_snp_counter/$snp_counter) * 100) . "%)") . "\n";

%seen = ();
@snps_gmap_alignment = grep { ! $seen{$_}++ } @snps_gmap_alignment;
my $snps_gmap_alignment_counter = 0;
foreach my $query_id (sort @snps_gmap_alignment){
	foreach my $snps_gmap_alignment_entry (@{$snps_positions->{$query_id}}){
		$snps_gmap_alignment_counter++;
	}
}

%seen = ();
@snps_no_gmap_alignment = grep { ! $seen{$_}++ } @snps_no_gmap_alignment;
my @snps_no_gmap_align = ();
my $snps_no_gmap_alignment_counter = 0;
foreach my $query_id (sort @snps_no_gmap_alignment){
	foreach my $snps_no_gmap_alignment_entry (@{$snps_positions->{$query_id}}){
		$snps_no_gmap_alignment_counter++;
		push(@snps_no_gmap_align, join("\t", $query_id, $snps_no_gmap_alignment_entry));
	}
}

my $total_gmap_snps_counter = ($snps_gmap_alignment_counter + $snps_no_gmap_alignment_counter);
print OUTFILE "Total SNP positions: $total_gmap_snps_counter" . "\n";
print OUTFILE "GMAP SNP positions: $snps_gmap_alignment_counter" . "\n";
print OUTFILE "SNP positions without a GMAP alignment: $snps_no_gmap_alignment_counter" . "\n";

my $snps_no_gmap_align_outfile = join('/', $output_dir, join("_", $snp_filename, "snps_no_gmap_align_outfile.txt"));
open(OUTFILE1, ">$snps_no_gmap_align_outfile") or die "Couldn't open file $snps_no_gmap_align_outfile for writting, $!";
foreach my $snps_no_gmap_align_entry (sort @snps_no_gmap_align){
	print OUTFILE1 $snps_no_gmap_align_entry . "\n";
}
close(OUTFILE1) or die "Couldn't close file $snps_no_gmap_align_outfile";

%seen = ();
@gmap_no_snps = grep { ! $seen{$_}++ } @gmap_no_snps;
my $gmap_no_snps_counter = scalar(@gmap_no_snps);

print OUTFILE "GMAP alignments without SNP positions: $gmap_no_snps_counter" . "\n";

my $gmap_no_snps_outfile = join('/', $output_dir, join("_", $snp_filename, "gmap_no_snps_outfile.txt"));
open(OUTFILE2, ">$gmap_no_snps_outfile") or die "Couldn't open file $gmap_no_snps_outfile for writting, $!";
foreach my $query_id (sort @gmap_no_snps){
	print OUTFILE2 $query_id . "\n";
}
close(OUTFILE2) or die "Couldn't close file $gmap_no_snps_outfile";
close(OUTFILE) or die "Couldn't close file $snps_stats_outfile";

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
				$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,$query_gap_position_list,
				$mismatch_position_list) = @split_gmap_align_entry;

	#  		  warn "$alignment_name\n";

			my $gmap_align_entry = join("\t", $coverage,$percent_identity,$matches,$mismatches,$indels,$unknowns,$query_id,
				$query_start,$query_end,$query_length,$query_strand,$target_id,$target_start,$target_end,$target_length,
				$target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,
				$query_align_block,$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,
				$query_gap_position_list,$mismatch_position_list);

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

	my (%fasta_seqs, %fasta_lengths);
	my $seqio = Bio::SeqIO->new(-file => $fasta_infile, '-format' => 'Fasta');
	while(my $seq_entry = $seqio->next_seq) {

		my $seq_id = $seq_entry->id;
		my $sequence = $seq_entry->seq;

	# 	    warn $seq_id . "\n";
	# 	    warn $sequence . "\n";
		$fasta_seqs{$seq_id} = $sequence;
		$fasta_lengths{$seq_id} = length($sequence);
	}

	# Clean out the sequence I/O object.
	$seqio = ();

	return (\%fasta_seqs, \%fasta_lengths);
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

	my $query_strand = shift;
	die "Error lost input query strand" unless defined $query_strand;

	my $amino_acid_start = shift;
	die "Error lost input amino acid start position" unless defined $amino_acid_start;

	my $snp_position = shift;
	die "Error lost input snp position" unless defined $snp_position;

	my $snp_variants = shift;
	die "Error lost input snp variants" unless defined $snp_variants;

	my ($snp_variant1, $snp_variant2) = split(/\//, $snp_variants);
	#       warn "BEFORE minus strand change " . join("/", $snp_variant1, $snp_variant2) . "\n";
	($snp_variant1, $snp_variant2) = (get_comp_base($snp_variant1), get_comp_base($snp_variant2)) if($query_strand eq "-");
	#       warn "AFTER minus strand change " . join("/", $snp_variant1, $snp_variant2) . "\n";

	my $adj_snp_position = $snp_position - $amino_acid_start;

	my ($snp_variants_list, $codon_list, $amino_acid_list, $substitution_type, $aa_sequence);

	my $subsitution_position = 0;
	my $is_snp_found = "false";

	my $amino_acid_position = 0;
	for(my $i = 0; $i < (length($dna_sequence) - 2); $i += 3) {
		if(($i > $adj_snp_position) and ($is_snp_found eq "false")){
			my $snp_codon_position = ($snp_position - (($amino_acid_start + $i) - 3));
			
			$snp_codon_position = abs($snp_codon_position - 2) if($query_strand eq "-");
			warn join("\t", $snp_position,(($amino_acid_start + $i) - 3), $snp_codon_position) . "\n";

			my $snp_codon = substr($dna_sequence, ($i - 3), 3);
	# 		  warn "BEFORE minus strand change " . $snp_codon . "\n";

			$snp_codon = reverse_complement($snp_codon) if($query_strand eq "-");
	# 		  warn "AFTER minus strand change " . $snp_codon . "\n";

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

	# 	    warn "BEFORE minus strand change " . $codon . "\n";

		$codon = reverse_complement($codon) if($query_strand eq "-");
	# 	    warn "AFTER minus strand change " . $codon . "\n";

		$aa_sequence .= convert_codon($codon);
		$amino_acid_position++;
	}
	my @split_aa_sequence = split('', $aa_sequence);
	my $amino_acid_base = $split_aa_sequence[$subsitution_position];
	$split_aa_sequence[$subsitution_position] = "[$amino_acid_base]";
	$split_aa_sequence[$subsitution_position] = "]" . $amino_acid_base. "[" if($query_strand eq "-");

	my $amino_acid_sequence = join('', @split_aa_sequence);
	$amino_acid_sequence = reverse($amino_acid_sequence) if($query_strand eq "-");
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

sub complement{

	my $dna_sequence = shift;
	die "Error lost input dna sequence" unless defined $dna_sequence;

	my @comp_seq = ();
	my @dna_seq = split('', $dna_sequence);
	for(my $i = 0; $i < length($dna_sequence); $i++) {
		my $nuc_base = $dna_seq[$i];
		my $comp_base = get_comp_base($nuc_base);
		push(@comp_seq, $comp_base);
	}


	return @comp_seq;
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
