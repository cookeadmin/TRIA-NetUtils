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
my ($input_dir, $output_dir, $infile_suffix, $hclust_suffix, $distFunction, $hclustFunction);
my @options = (
	'i=s',	\$input_dir,
	'o=s',	\$output_dir,
	's=s',	\$infile_suffix,
	't=s',	\$hclust_suffix,
	'd=s',	\$distFunction,
	'h=s',	\$hclustFunction
);
&GetOptions(@options);

usage() unless (
	defined $input_dir
	and defined $output_dir
);
# Parameter options initialized if not defined in @options argument strings.

# R input file suffix. Can be anything (.tsv, etc.)
$infile_suffix = 'tsv' unless defined $infile_suffix;

# hclust input file suffix. Can be anything (.hclust, etc.)
$hclust_suffix = 'hclust' unless defined $hclust_suffix;

# The distance measure to be used on the data. This should be (an unambiguous abbreviation of or full word) one of either "euclidean", "maximum", "manhattan" or "binary". The default is "euclidean".
$distFunction = 'euclidean' unless defined $distFunction;

# The agglomeration method to be used in hierarchial clustering. This should be (an unambiguous abbreviation of or full word) one of either "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". The default is "average".
$hclustFunction = 'average' unless defined $hclustFunction;


# Usage details on input parameters for the hierarchial clustering of gene ids between multiple samples.
sub usage {

die <<"USAGE";

Usage: $0 -i input_dir -o output_dir -s infile_suffix -t hclust_suffix -d distFunction -h hclustFunction

OPTIONS:

-i input_dir - The input_dir containing files in the following format. Where the header consists of the sample names to compare and each subsequent line indicates the gene ids and the log2 fold change values found within the samples specified in the header line.
e.g.
ID_REF	JP.1dpi.WDF.vs.WDC	JP.1dpi.WWF.vs.WWC	JP.28dpi.WDC.vs.WWC	JP.28dpi.WDF.vs.WDC	JP.28dpi.WDF.vs.WDW	JP.28dpi.WDF.vs.WWF	JP.28dpi.WDW.vs.WWW	JP.28dpi.WWF.vs.WWC	JP.28dpi.WWF.vs.WWW	JP.7dpi.WDC.vs.WWC	JP.7dpi.WDF.vs.WDC	JP.7dpi.WDF.vs.WDW	JP.7dpi.WDF.vs.WWF	JP.7dpi.WWF.vs.WWC	JP.7dpi.WWF.vs.WWW	JP.vs.LP.1dpi.WDF	JP.vs.LP.1dpi.WWF	JP.vs.LP.28dpi.WDF	JP.vs.LP.28dpi.WWF
1.1.13	0	0.731183242	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0

-o output_dir - The directory to output the hclust associated images and files.

-s infile_suffix - input file suffix. Can be anything (.tsv, etc.)

-t hclust_suffix - hclust input file suffix. Can be anything (.hclust, etc.)

-d distFunction - The distance measure to be used on the data. This should be (an unambiguous abbreviation of or full word) one of either "euclidean", "maximum", "manhattan" or "binary". The default is "euclidean".

-h hclustFunction - The agglomeration method to be used in hierarchial clustering. This should be (an unambiguous abbreviation of or full word) one of either "ward", "single", "complete", "average", "mcquitty", "median" or "centroid". The default is "average".

USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

# Parse the R input file directory for any R input files with the specified suffix.
my $input_files = find_input_files($input_dir,$infile_suffix);

# Foreach input file in the R input file list perform hierarchial clustering on gene ids between multiple samples.
foreach my $input_file (@{$input_files}){
#       warn "$input_file\n";
      # Get the R input file path.
      my $infile = join('/', $input_dir, $input_file);

      # Extract sample names from the input filename for the hclust dendrogram plot title and the outfile name.
      my ($title, $outfile);
      if($input_file =~ m/([\w\_\-\.]+)\.$infile_suffix/){
	    $outfile = $1;

	    $title = $outfile;
      }

      my $r_infile = join('/', $output_dir, $input_file);
      my $i = 0;
      open INFILE, "<$infile" or die "Error opening $infile for reading: $!";
      open OUTFILE, ">$r_infile" or die "Error opening $r_infile for writting: $!";
      while(<INFILE>){
	    chomp $_;
	    my @split_entry = split(/\t/, $_);
	    if($i eq 0){
		  
		  print OUTFILE join("\t", @split_entry[1..$#split_entry]) . "\n";
	    }
	    else{
		  print OUTFILE $_ . "\n";
	    }
	    $i++;
	    

      }
      close INFILE or die "Error closing $infile: $!";
      close OUTFILE or die "Error closing $r_infile: $!";

      # The hclust sample names output filename.
      my $hclust_samples_outfile = join('/', $output_dir, $outfile . "-sample-names");

      # The hclust gene ids output filename.
      my $hclust_gene_ids_outfile = join('/', $output_dir, $outfile . "-gene-ids");

# Pipe the following R script that computes hierarchial clustering method on the R input file to the R program.
open (R_SCRIPT_PIPE, "|R --vanilla --slave");
print R_SCRIPT_PIPE <<EOF;
options(warn=1)
## Loading Global library packages.

## Cluster and Tree Conversion Package for writting Newick Tree files.
library(ctc)

## Read the data file into a data frame.
dataframe <- read.delim("$r_infile")

datacolumns <- names(dataframe)
datamatrix <- as.matrix(dataframe[datacolumns])

## Sample names are the transposed matrix
samples <- t(datamatrix)

## Hierarchical Clustering
## Generate the sample names hclust dendrogram.
hclustSamples <- hclust(dist(samples, method = "$distFunction"), method="$hclustFunction")
png('$hclust_samples_outfile-hclust-dendrogram.png\');
plot(hclustSamples, hang = -1, main="$title")

## Generate sample names table from hclust results.
hclustSamplesDendroLabels <- data.frame(labels=hclustSamples\$labels[hclustSamples\$order])
write.table(hclustSamplesDendroLabels, file = "$hclust_samples_outfile.$hclust_suffix", sep = "\t")

## Generate Samples newick tree from HCLUST results.
write(hc2Newick(hclustSamples),file='$hclust_samples_outfile-hclust.newick')

## Turn off R graphical device
## If successful, prints on STDOUT
## null device
##          1
dev.off()

## Gene ids are the direct matrix
geneIds <- (datamatrix)

## Hierarchical Clustering
## Generate the gene ids hclust dendrogram.
hclustGeneIds <- hclust(dist(geneIds, method = "$distFunction"), method="$hclustFunction")
# png('$hclust_gene_ids_outfile-hclust-dendrogram.png\');
# plot(hclustGeneIds, hang = -1, main="$title")

## Generate gene ids table from hclust results.
hclustGeneIdsDendroLabels <- data.frame(labels=hclustGeneIds\$labels[hclustGeneIds\$order])
write.table(hclustGeneIdsDendroLabels, file = "$hclust_gene_ids_outfile.$hclust_suffix", sep = "\t")

## Generate OTUs newick tree from HCLUST results.
write(hc2Newick(hclustGeneIds),file='$hclust_gene_ids_outfile-hclust.newick')

## Turn off R graphical device
## If successful, prints on STDOUT
## null device
##          1
#dev.off()

## Exit out of R program
q()

EOF
close(R_SCRIPT_PIPE);

      # Parse the hclust input file directory for any hclust input files with the specified suffix.
      my $hclust_input_files = find_input_files($output_dir, $hclust_suffix);

      # Foreach input file in the hclust input file list, parse the hclust gene id or sample name lists for a cleaner file.
      foreach my $filename (@{$hclust_input_files}){

	    # Get the hclust input file path.
	    my $hclust_infile = join('/', $output_dir, $filename);

	    my $outfile;
	    if($filename =~ m/([\w\_\-\.]+)\.$hclust_suffix/){
		  $outfile = $1;
	    }
	    my $hclust_outfile = join('/', $output_dir, $outfile . "-hclust-dendrogram.txt");
	    open INFILE, "<$hclust_infile" or die "Error opening $hclust_infile for reading: $!";
	    open OUTFILE, ">$hclust_outfile" or die "Error opening $hclust_outfile for writting: $!";
	    while(<INFILE>){
		  chomp $_;
		  if ($_ =~ m/^"\d+"\t"(.+)"$/){
			      my $clustered_item = $1;
      # 			$clustered_item =~ s/\./-/g;
			      print OUTFILE "$clustered_item\n";
		  }
		  
	    }
	    close INFILE or die "Error closing $hclust_infile: $!";
	    close OUTFILE or die "Error closing $hclust_outfile: $!";
	    if(-s $hclust_outfile){
		  unlink $hclust_infile;
	    }
      }

      my @hclust_gene_ids;
      my $hclust_gene_ids_infile = $hclust_gene_ids_outfile . "-hclust-dendrogram.txt";
      open INFILE, "<$hclust_gene_ids_infile" or die "Error opening $hclust_gene_ids_infile for reading: $!";
      while(<INFILE>){
	    chomp $_;
      #       warn "$_\n";
	    push(@hclust_gene_ids, $_);

      }
      close INFILE or die "Error closing $hclust_gene_ids_infile: $!";

      my @hclust_sample_names;
      my $hclust_samples_infile = $hclust_samples_outfile . "-hclust-dendrogram.txt";
      open INFILE, "<$hclust_samples_infile" or die "Error opening $hclust_samples_infile for reading: $!";
      while(<INFILE>){
	    chomp $_;
      #       warn "$_\n";
	    push(@hclust_sample_names, $_);

      }
      close INFILE or die "Error closing $hclust_samples_infile: $!";

      open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
      $i = 0;
      my (%column_index, %column_values);
      while(<INFILE>){
	    chomp $_;
	    my @split_entry = split(/\t/, $_);
	    if($i eq 0){
		  for(my $j = 1; $j < scalar(@split_entry); $j++){
			$column_index{$j} = $split_entry[$j];
		  }
	    }
	    else{
		  for(my $j = 1; $j < scalar(@split_entry); $j++){
			my ($gene_ids, $sample_names) = ($split_entry[0], $column_index{$j});
			$column_values{$gene_ids}{$sample_names} = $split_entry[$j];
		  }
	    }
	    $i++;
      }
      close(INFILE) or die "Couldn't close file $infile";

      my $r_outfile = join('/', $output_dir, $outfile . "-hclust.tsv");
      open OUTFILE, ">$r_outfile" or die "Error opening $r_outfile for writting: $!";
      print OUTFILE join("\t", @hclust_sample_names) . "\n";
      my $gene_id_counter = 0;
      my @log_FC_values;
      foreach my $gene_ids (@hclust_gene_ids){
	    my @row_entry;
	    foreach my $sample_names (@hclust_sample_names){
		  push(@row_entry, $column_values{$gene_ids}{$sample_names});
		  push(@log_FC_values, $column_values{$gene_ids}{$sample_names});
	    }
	    print OUTFILE join("\t", $gene_ids, @row_entry) . "\n";
	    $gene_id_counter++;
      }
      close OUTFILE or die "Error closing $r_outfile: $!";

      my $min_FC_value = floor(min(@log_FC_values));
      my $max_FC_value = ceil(max(@log_FC_values));

      # The heatmap output image filename.
      my $heatmap_overview_file = join('/', $output_dir, $outfile . "-heatmap-overview");
      # Pipe the following R script to the R program that generates a heatmap using the R input file.
      # open (R_SCRIPT_PIPE, "|R --vanilla --slave");
open (R_SCRIPT_PIPE, "|R --vanilla --slave");
print R_SCRIPT_PIPE <<EOF;
options(warn=1)
## R Color Brewer Package
library(RColorBrewer)
## Pretty Heatmap Package
library(pheatmap)

## Read the data file into a data frame.
datamatrix <- as.matrix(read.delim("$r_infile", header = TRUE,sep="\t"))
dim(datamatrix)

## Generate Heatmap
pheatmap(datamatrix,
filename = "$heatmap_overview_file.pdf",
border_color = "#202020",
cellwidth = 60,
cellheight = 12,
color = colorRampPalette(c("red","black","green"))(299),
cluster_rows = TRUE,
cluster_cols = TRUE,
clustering_distance_cols = "$distFunction",
clustering_distance_rows = "$distFunction",
clustering_method = "$hclustFunction",
bg = "white",
treeheight_row = 0,
treeheight_col = 500,
legend = TRUE,
breaks = unique(c( seq($min_FC_value,-0.1,length=100),seq(-0.1,0.1,length=100),seq(0.1,$max_FC_value,length=100))),
compress = TRUE,
)

## Turn off R graphical device
## If successful, prints on # STDOUT
dev.off()

## Exit out of R program
q()
EOF
close(R_SCRIPT_PIPE);
	


  
      # Create output directory if it doesn't already exist.
      my $heatmap_output_dir = join('/', $output_dir, "heatmap_images");
      unless(-d $heatmap_output_dir){
	    mkdir($heatmap_output_dir, 0777) or die "Can't make directory: $!";
      }
      
      system("convert -density 100 $heatmap_overview_file.pdf $heatmap_overview_file.png"
      ) == 0 or die "Error calling convert: $?";

      my $heatmap_image = Image::Magick->new;
      $heatmap_image->Read("$heatmap_overview_file.png");
      
      my $crop_geometry = join("", $heatmap_image->Get('width'), "x675+0+0");

      my $heatmap_dendrogram = join('/', $heatmap_output_dir, $outfile . "-heatmap-dendrogram.png");

      system("convert -density 100 $heatmap_overview_file.png -crop $crop_geometry $heatmap_dendrogram"
      ) == 0 or die "Error calling convert: $?";

      my $gene_id_remainder = ($gene_id_counter % 1000);
      my $gene_id_quotient = ($gene_id_counter - $gene_id_remainder);

      my $files_expected = $gene_id_quotient / 1000;
# 	    if(($i < $gene_id_quotient) and ($file_counter < $files_expected)){
# 		  ($start_index, $end_index) = (($i+1), ($i+1000));
# 	    }elsif(($i < $gene_id_quotient) and ($gene_id_remainder <= 100) and ($file_counter eq $files_expected)){
# 		  ($start_index, $end_index) = (($i+1), ($i+1000+$gene_id_remainder));
# 	    }elsif(($i eq $gene_id_quotient) and ($gene_id_remainder > 100) and ($file_counter eq $files_expected)){
# 		  ($start_index, $end_index) = (($i+1), ($i+$gene_id_remainder));
# 	    }

      
      my $file_counter = 1;
      for (my $i = 0; $i <= $gene_id_quotient; $i += 1000){

	    my ($start_index, $end_index);
	    if($i < $gene_id_quotient){
		  if(($gene_id_remainder <= 100) and ($file_counter eq $files_expected)){
			($start_index, $end_index) = (($i+1), ($i+1000+$gene_id_remainder));
		  }else{
			($start_index, $end_index) = (($i+1), ($i+1000));
		  }
	    }elsif(($i eq $gene_id_quotient) and ($gene_id_remainder > 100) and ($file_counter > $files_expected)){
		  ($start_index, $end_index) = (($i+1), ($i+$gene_id_remainder));
	    }else{
		last; # we want to break out of the loop because we are done.
	    }

	    # The heatmap output image filename.
	    my $heatmap_outfile = join('/', $heatmap_output_dir, $outfile . "-heatmap-$start_index-$end_index-ids");

	    # Pipe the following R script to the R program that generates a heatmap using the R input file.
	    # open (R_SCRIPT_PIPE, "|R --vanilla --slave");
warn "Generating $outfile-heatmap-$start_index-$end_index-ids.png\n";
open (R_SCRIPT_PIPE, "|R --vanilla --slave");
print R_SCRIPT_PIPE <<EOF;
options(warn=1)
## R Color Brewer Package
library(RColorBrewer)
## Pretty Heatmap Package
library(pheatmap)

## Read the data file into a data frame.
datamatrix <- as.matrix(read.delim("$r_outfile", header = TRUE,sep="\t"))
dim(datamatrix)
## Generate Heatmap
pheatmap(datamatrix[$start_index:$end_index, ],
filename = "$heatmap_outfile.pdf",
border_color = "#202020",
cellwidth = 60,
cellheight = 12,
color = colorRampPalette(c("red","black","green"))(299),
cluster_rows = FALSE,
cluster_cols = FALSE,
# clustering_distance_cols = "$distFunction",
# clustering_distance_rows = "$distFunction",
# clustering_method = "$hclustFunction",
bg = "white",
treeheight_row = 0,
treeheight_col = 500,
legend = TRUE,
breaks = unique(c( seq($min_FC_value,-0.1,length=100),seq(-0.1,0.1,length=100),seq(0.1,$max_FC_value,length=100))),
compress = TRUE,
)

## Turn off R graphical device
## If successful, prints on STDOUT
dev.off()
	

## Exit out of R program
q()
EOF
close(R_SCRIPT_PIPE);

	    $file_counter++;
	    system("convert -density 100 $heatmap_outfile.pdf $heatmap_outfile.png"
	    ) == 0 or die "Error calling convert: $?";

	    system("convert $heatmap_dendrogram $heatmap_outfile.png -append $heatmap_outfile.png"
	    ) == 0 or die "Error calling convert: $?";
	    unlink "$heatmap_outfile.pdf";
      }
      unlink "$heatmap_overview_file.png";
      unlink $heatmap_dendrogram;
}



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
