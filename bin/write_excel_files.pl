#!/usr/bin/perl
use warnings;
use strict;
use Spreadsheet::WriteExcel;

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

sub usage {
    
die << "USAGE";

    Usage: $0 -i input_dir -o output_dir

USAGE
}
## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

write_excel_files($input_dir, $output_dir);

sub write_excel_files{
    
      my $input_dir = shift;
      die "Error lost input directory" unless defined $input_dir;

      my $output_dir = shift;
      die "Error lost output directory" unless defined $output_dir;


      my $files = find_files($input_dir, "txt");

      foreach my $file (keys %{$files}){
	    my $filename = $files->{$file};

	    my @excel_output;
	    open(INFILE, "<$filename") or die "Couldn't open file $filename for reading, $!";
	    while(<INFILE>){
		  chomp $_;
		  #warn $_ . "\n";
		  push(@excel_output, $_);

	    }
	    close(INFILE) or die "Couldn't close file $filename";

	    my $excel_filename = $file;

	    $excel_filename =~ s/txt/xls/g;

	    warn "Writting excel file $excel_filename....." . "\n";
	    my $excel_file = join('/', $output_dir, $excel_filename);

	    my $workbook  = Spreadsheet::WriteExcel->new($excel_file);
	    my $worksheet = $workbook->add_worksheet();
	    my $i = 0;
	    foreach my $row_entry (@excel_output){
		  
		  my @split_row = split(/\t/, $row_entry);
		  my $j = 0;
		  foreach my $entry (@split_row){
			$worksheet->write($i, $j,  $entry);
			$j++;
		  }
		  $i++;
	    }

      }

}

sub find_files {
	my $dir = shift;
	die "Error lost directory" unless defined $dir;
	my $suffix = shift;
	die "Error lost pattern" unless defined $suffix;
    
    #die $dir ."\n";
	opendir(DIR,$dir) or die "Error opening $dir: $!";

	my %files;
	while(my $file = readdir(DIR)){
	    my $filename = join('/', $dir, $file);
	    $files{$file} = $filename if ($file =~ /\.$suffix$/);
	}
	closedir(DIR) or die "Error closing $dir: $!";
	return \%files;
}