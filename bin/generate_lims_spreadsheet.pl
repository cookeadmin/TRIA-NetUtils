#!/usr/bin/perl
use warnings;
use strict;
use Spreadsheet::WriteExcel;

use Getopt::Long;

my $output_dir;

my @options = (
    'o=s'    => \$output_dir,
);
&GetOptions(@options);

usage() unless (
    defined $output_dir
);

sub usage {
    
die << "USAGE";

    Usage: $0 -o output_dir

USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}



# my $excel_filename = $file;
# 
# 
# warn "Writting excel file $excel_filename....." . "\n";
# my $excel_file = join('/', $output_dir, $excel_filename);
# 
# my $workbook  = Spreadsheet::WriteExcel->new($excel_file);
# my $worksheet = $workbook->add_worksheet();
# my $i = 0;
# foreach my $row_entry (@excel_output){
#       
#       my @split_row = split(/\t/, $row_entry);
#       my $j = 0;
#       foreach my $entry (@split_row){
# 	    $worksheet->write($i, $j,  $entry);
# 	    $j++;
#       }
#       $i++;
# }


sub get_lims_sheets{
      my %lims_sheets;
      my @lims_sheet_names = ("Experiment Description", "Experiment Factors", "Biomaterials");
      @{$lims_sheets{"Experiment Description"}} = ("*Project", "*Contact person", "*Code", "Category", "Name", "Description", "Comment", "Start date", "Location");
      @{$lims_sheets{"Experiment Factors"}} = ("Factor", "*Name", "Description", "Variables");
      @{$lims_sheets{"Biomaterials"}} = ("*Species", "*Clone/Family", "*Organism description", "*Details of culture or collection", "", "");

}

