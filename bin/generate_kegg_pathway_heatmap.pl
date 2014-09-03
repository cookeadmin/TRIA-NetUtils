#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Image::Magick;
use List::Util qw( min max );
use POSIX qw( floor ceil );
# Program for the hierarchial clustering of gene ids between multiple samples.

# Command line example for supplying program with the proper parameter options.
# perl generate_heatmap.pl -i /path/to/heatmap_input_R/JP-DE-all-days-comparison/submission-summary -o /path/to/heatmap_output_R/JP-DE-all-days-comparison

# Declaration of parameter options to be initialized in @options argument strings.
my ($input_image, $output_dir);
my @options = (
	'i=s',	\$input_image,
	'o=s',	\$output_dir
);
&GetOptions(@options);

# usage() unless (
# # 	defined $input_image
# # 	and defined $output_dir
# );
# Parameter options initialized if not defined in @options argument strings.


# Usage details on input parameters for the hierarchial clustering of gene ids between multiple samples.
sub usage {

die <<"USAGE";

Usage: $0 -i input_dir -o output_dir -s infile_suffix -t hclust_suffix -d distFunction -h hclustFunction

OPTIONS:

-i input_image - The input_dir containing files in the following format. Where the header consists of the sample names to compare and each subsequent line indicates the gene ids and the log2 fold change values found within the samples specified in the header line.

-o output_dir - The directory to output the hclust associated images and files.

USAGE
}

# Create output directory if it doesn't already exist.
# unless(-d $output_dir){
# 	mkdir($output_dir, 0777) or die "Can't make directory: $!";
# }

open (R_SCRIPT_PIPE, "|R --vanilla --slave");
print R_SCRIPT_PIPE <<EOF;
library(png)
img <- readPNG("/home/cookeadmin/workspace/adriana/kegg_pathway/map00902.png")
image(as.raster(img))

## Turn off R graphical device
## If successful, prints on STDOUT
## null device
##          1
dev.off()

## Exit out of R program
q()

EOF
close(R_SCRIPT_PIPE);


###
exit 0;

=head2 B<\@input_files> = B<find_input_files(B<$input_dir>, B<$suffix>)> - Find all input files with the directory and with the suffix specified and returns a reference to the list of input files

I<Input Parameters>:

B<$input_dir>  - The input file directory to parse the input files.

B<$suffix>  - Input file suffix. Can be anything (.tsv, etc.)

I<Output Parameters>:

B<\@input_files> - Reference to a list of input files.

E<10>

E<10>

=cut
sub find_input_files {

	# The directory to parse for input files.
	my $input_dir = shift;
	die "Error lost input file directory" unless defined $input_dir;

	# Input file suffix. Can be anything (.tsv, etc.)
	my $suffix = shift;
	die "Error lost input file suffix" unless defined $suffix;

	if(-d $input_dir){ # Check if $input_dir is a directory.
		opendir(DIR,$input_dir) or die "Error opening $input_dir: $!";
		
		warn "Loading the following files.....\n";

		my @input_files;
		while(my $input_file = readdir(DIR)){
			if ($input_file =~ /\.$suffix$/){
			      warn $input_file . "\n";
			      push(@input_files, $input_file);
			}
		}
		closedir(DIR) or die "Error closing $input_dir: $!";

		# Reference to a list of input files.
		return \@input_files;
	}
	else{
		die "Error $input_dir does not exist!\n";
	}
}
