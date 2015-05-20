#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my ($input_dir, $gff3_filename, $min_percent_id, $output_dir);
GetOptions(
      'i=s'    => \$input_dir,
      'n=s'    => \$gff3_filename,
      'p=s'    => \$min_percent_id,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $input_dir
      and $gff3_filename
      and defined $output_dir
);

# The fastq input file type. Default: gzfastq
$min_percent_id = 80 unless defined $min_percent_id;

sub usage {

die <<"USAGE";

Usage: $0 -i input_dir -n gff3_filname -p min_percent_id -o output_dir

Description - 

OPTIONS:

      -i input_dir -
      -n gff3_filname - 
      -p min_percent_id
      -o output_dir -

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my ($blat2gff3_files, $blat2gff3_file_count) = find_files($input_dir, "gff");

my %blat2gff3_hits = ();
foreach my $file_name (sort keys %{$blat2gff3_files}){
	warn "Processing " . $file_name . ".....\n";
	my $blat2gff3_infile = $blat2gff3_files->{$file_name};
	
	my $i = 0;
	open(INFILE,"<$blat2gff3_infile") or die "Error opening $blat2gff3_infile: $!";
	while(<INFILE>){
		chomp $_;
		if($i ne 0){
			warn $_ . "\n";
			my @split_blastx_entry = split(/\t/, $_);
			
			# 142155638       BLAT    cDNA_match      182     2326    98.48   -       .       ID=GQ03808_D10.1_mid2;Target=GQ03808_D10.1 1121 1
			my ($target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_blastx_entry;
			if($match_type eq "cDNA_match"){
				my $blat2gff3_entry = join("\t", $target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
				
				my $query_id = "";
				if($attributes =~ m/ID=(.+);Target=.+\s\d+\s\d+/){ # ID=WS00720_I24.1_mid7;Target=WS00720_I24.1 908 30
					$query_id = $1;
				}else{
					die "Error: Attributes not in correct format: $attributes";
				}
				
				my $unique_id = join("_", $target_id, $query_id);
				$blat2gff3_hits{$unique_id}{"cDNA_match"} = $blat2gff3_entry;
			}
			if($match_type eq "exon"){
				my $blat2gff3_entry = join("\t", $target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
				
				my $query_id = "";
				if($attributes =~ m/Parent=(.+);Target=.+\s\d+\s\d+/){ # Parent=WS0337_D11.1_mid1;Target=WS0337_D11.1 20 6
					$query_id = $1;
				}else{
					die "Error: Attributes not in correct format: $attributes";
				}
				
				my $unique_id = join("_", $target_id, $query_id);
				push(@{$blat2gff3_hits{$unique_id}{"exon"}}, $blat2gff3_entry);
			}
		}
		$i++;
	}
	close INFILE or die "Error closing $blat2gff3_infile: $!";
}

my %filtered_blat2gff3 = ();
foreach my $unique_id (sort {$a cmp $b} keys %blat2gff3_hits){
	my @split_blastx_entry = split(/\t/, $blat2gff3_hits{$unique_id}{"cDNA_match"});
	my ($target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_blastx_entry;
	if($percent_id >= $min_percent_id){
		my $query_id;
		if($attributes =~ m/ID=.+;Target=(.+)\s\d+\s\d+/){ # ID=WS00720_I24.1_mid7;Target=WS00720_I24.1 908 30
			$query_id = $1;
		}else{
			die "Error: Attributes not in correct format: $attributes";
		}
		my @blat2gff3_entries = ();
		push(@blat2gff3_entries, $blat2gff3_hits{$unique_id}{"cDNA_match"});
		my $i = 1;
		foreach my $exon_entry (@{$blat2gff3_hits{$unique_id}{"exon"}}){
		
			my @split_blastx_entry = split(/\t/, $exon_entry);
			my ($target_id, $program_source, $match_type, $target_start, $target_end, $percent_id, $strand, $frame, $attributes) = @split_blastx_entry;
			my $exon_num = join("_", $match_type, $i);
			
			my $new_exon_entry = join("\t", $target_id, $program_source, $exon_num, $target_start, $target_end, $percent_id, $strand, $frame, $attributes);
			push(@blat2gff3_entries, $new_exon_entry);
			$i++;
		}
		
		my $new_blat2gff3_entry = join("\n", @blat2gff3_entries);
		
 		push(@{$filtered_blat2gff3{$target_id}{$percent_id}}, $new_blat2gff3_entry);
	}
}

my $blat2gff_outfile = join('/', $output_dir, join("-", $gff3_filename, "query-id-sorted.gff"));
open(OUTFILE, ">$blat2gff_outfile") or die "Couldn't open file $blat2gff_outfile for writting, $!";
print OUTFILE "##gff-version 3\n";
foreach my $query_id (sort {$a cmp $b} keys %filtered_blat2gff3){
	foreach my $percent_id (sort {$b <=> $a} keys $filtered_blat2gff3{$query_id}){
		foreach my $blat2gff3_entry (@{$filtered_blat2gff3{$query_id}{$percent_id}}){
			print OUTFILE $blat2gff3_entry . "\n";
		}
	}
}
close(OUTFILE) or die "Couldn't close file $blat2gff_outfile";

my $blat2gff_outfile = join('/', $output_dir, join("-", $gff3_filename, "contig-sorted.gff"));
open(OUTFILE, ">$blat2gff_outfile") or die "Couldn't open file $blat2gff_outfile for writting, $!";
print OUTFILE "##gff-version 3\n";
foreach my $target_id (sort {$a cmp $b} keys %filtered_blat2gff3){
	foreach my $percent_id (sort {$b <=> $a} keys $filtered_blat2gff3{$target_id}){
		foreach my $blat2gff3_entry (@{$filtered_blat2gff3{$target_id}{$percent_id}}){
			print OUTFILE $blat2gff3_entry . "\n";
		}
	}
}
close(OUTFILE) or die "Couldn't close file $blat2gff_outfile";
# (\%files, $file_counter) = find_files($infile_dir) - Find all files in the specified input file directory with the file extension *.suffix.
#
# Input paramater(s):
#
# $infile_dir - The input file directory.
#
# $suffix - The file extension suffix.
#
# Output paramater(s):
#
# \%files - A hash reference containing all the files with file extension *.suffix in key/value pairs.
#
# key => filename ( e.g. filename.suffix )
# value => absolue filepath ( e.g. /path/to/filename.suffix )
#
# $file_count - The number of files stored with file extension *.suffix.
sub find_files{
    
	# The input file directory.
	my $infile_dir = shift;
	die "Error lost input file directory" unless defined $infile_dir;
	
	# The file extension suffix.
	my $suffix = shift;
	die "Error lost file extension suffix directory" unless defined $suffix;
	
	if(-d $infile_dir){ # Check if $infile_dir is a directory.
		my %files = ();
		my $file_counter = 0;
		opendir(DIR, $infile_dir) || die "Error in opening dir $infile_dir\n";
		while( my $file_name = readdir(DIR)){
			my $infile_name = join('/', $infile_dir, $file_name) if ($file_name =~ m/\.$suffix$/);
			warn "$infile_name\n" if ($file_name =~ m/\.$suffix$/);
			$files{$file_name} = $infile_name if ($file_name =~ m/\.$suffix$/);
			$file_counter++ if ($file_name =~ m/\.$suffix$/);
		}
		closedir(DIR);
		
		return (\%files, $file_counter);
	}else{
		die "Error $infile_dir does not exist!\n";
	}
}
