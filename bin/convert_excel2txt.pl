#!/usr/bin/perl
use warnings;
use strict;

use Text::Iconv;
use Spreadsheet::XLSX;
use Spreadsheet::ParseExcel;
use Spreadsheet::Read;

use Getopt::Long;

my ($input_dir, $output_dir);

my @options = (
    'i=s'    => \$input_dir,
    'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
    defined $input_dir
    and $output_dir
);

#perl ncbi-microarray-submission.pl -i /home/cookeadmin/workspace/Adriana/GEO/Phloem/JP-vs-LP-microarray-data-final/JP-DE/differentially-expressed-genes-xls -o /home/cookeadmin/workspace/Adriana/GEO/Phloem/JP-vs-LP-microarray-data-final/JP-DE

sub usage {
    
die << "USAGE";

    Usage: $0 -i input_dir -o output_dir

    e.g. perl $0 -i  -o

    Description -

    OPTIONS:
    -i input_dir -
    -o output_dir -

USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}


my $outfile_dir = parse_xlsx_files($input_dir, $output_dir);

sub find_files {
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	my $suffix = shift;
	die "Error lost pattern" unless defined $suffix;
    
	opendir(DIR,$dir) or die "Error opening $dir: $!";

	my %files;
	while(my $file = readdir(DIR)){
	    my $filename = join('/', $dir, $file);
	    $files{$file} = $filename if ($file =~ /\.$suffix$/);
	}
	closedir(DIR) or die "Error closing $dir: $!";
	return \%files;
}

sub parse_xlsx_files{
    
    my $input_dir = shift;
    die "Error lost input directory" unless defined $input_dir;

    my $output_dir = shift;
    die "Error lost output directory" unless defined $output_dir;
    
    ## Create output directory if it doesn't already exist.

    # Text::Iconv is not really required.
    # This can be any object with the convert method. Or nothing.
    my $converter = Text::Iconv ->new("utf-8", "windows-1251");
    
    my $xlsx_files = find_files($input_dir, "xlsx");
    
    foreach my $xlsx_file (sort {$a cmp $b} keys %{$xlsx_files}){
	my $filename = $xlsx_file;
	$filename =~ s/\.xlsx/\.txt/g;
	
      print join(" ", "Parsing xlsx file:", $xlsx_file . "\n");
      
      my $excel = Spreadsheet::XLSX->new($xlsx_files->{$xlsx_file}, $converter);
      
      foreach my $sheet (@{$excel->{Worksheet}}){

	    my @excel_output = ();

	    $sheet->{MaxRow} ||= $sheet->{MinRow};
	    
	    foreach my $row ($sheet->{MinRow}..$sheet->{MaxRow}){
	    
		  $sheet->{MaxCol} ||= $sheet->{MinCol};
		  
		  my @row_values;
		  foreach my $col ($sheet->{MinCol}..$sheet->{MaxCol}){
			
			my $cell = $sheet->{Cells}[$row][$col];
			
			if($cell){
			      push(@row_values, $cell->{Val});
			}
			
		  }
		  my $row_entry = join("\t", @row_values);
		  push(@excel_output, $row_entry);
	    }
		  
	    my $outfile = join('/', $output_dir, $sheet->{Name} . ".csv");
	    open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
	    foreach my $entry (@excel_output){
	    
		  print OUTFILE $entry . "\n";
	    }
	    close(OUTFILE) or die "Couldn't close file $outfile"; 
	    
      }

	
    }
    return $outfile_dir;
}