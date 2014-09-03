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

$outfile_dir = parse_xls_files($input_dir, $output_dir);

## Create output directory if it doesn't already exist.
my $summary_dir = join('/', $output_dir, "submission-summary");
unless(-d $summary_dir){
	mkdir($summary_dir, 0777) or die "Can't make directory: $!";
}
parse_txt_files($outfile_dir, $summary_dir);


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

sub parse_xlsx_files{
    
    my $input_dir = shift;
    die "Error lost input directory" unless defined $input_dir;

    my $output_dir = shift;
    die "Error lost output directory" unless defined $output_dir;
    
    ## Create output directory if it doesn't already exist.
    my $outfile_dir = join('/', $output_dir, "differentially-expressed-genes-txt");
    unless(-d $outfile_dir){
	mkdir($outfile_dir, 0777) or die "Can't make directory: $!";
    }
    
    # Text::Iconv is not really required.
    # This can be any object with the convert method. Or nothing.
    my $converter = Text::Iconv ->new("utf-8", "windows-1251");
    
    my $xlsx_files = find_files($input_dir, "xlsx");
    
    foreach my $xlsx_file (sort {$a cmp $b} keys %{$xlsx_files}){
	my $filename = $xlsx_file;
	$filename =~ s/vs/v/g;
	$filename =~ s/v/vs/g;
	$filename =~ s/\.xlsx/\.txt/g;
	
	my $outfile = join('/', $outfile_dir, $filename);
	unless(-s $outfile){
	    print join(" ", "Parsing xlsx file:", $xlsx_file . "\n");
	    
	    my $excel = Spreadsheet::XLSX->new($xlsx_files->{$xlsx_file}, $converter);
	    
	    my @excel_output;
	    foreach my $sheet (@{$excel->{Worksheet}}){
		
		if($sheet->{Name} eq "Sheet1"){
		    
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
		    
		}
	    }

	    open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
	    foreach my $entry (@excel_output){
		
		print OUTFILE $entry . "\n";
	    }
	    close(OUTFILE) or die "Couldn't close file $outfile";
	}
    }
    return $outfile_dir;
}

sub parse_xls_files{
    
	my $input_dir = shift;
	die "Error lost input directory" unless defined $input_dir;
    
	my $output_dir = shift;
	die "Error lost output directory" unless defined $output_dir;
    
    ## Create output directory if it doesn't already exist.
    my $outfile_dir = join('/', $output_dir, "differentially-expressed-genes-txt");
    unless(-d $outfile_dir){
	mkdir($outfile_dir, 0777) or die "Can't make directory: $!";
    }
    
    my $xls_files = find_files($input_dir, "xls");
    foreach my $xls_file (sort {$a cmp $b} keys %{$xls_files}){
	my $infile = $xls_files->{$xls_file};
	my $filename = $xls_file;
	$filename =~ s/vs/v/g;
	$filename =~ s/v/vs/g;
	$filename =~ s/\.xls/\.txt/g;

	my $outfile = join('/', $outfile_dir, $filename);
	unless(-s $outfile){
	    
	    print join(" ", "Parsing xls file:", $xls_file . "\n");
	    my $parser = Spreadsheet::ParseExcel->new();
	    my $workbook = $parser->parse($infile);
	    
	    if (!defined $workbook){
		die $parser->error(), ".\n";
	    }
	    
	    my @excel_output;
	    for my $worksheet ($workbook->worksheets()){
		
		if($worksheet->get_name() eq "Sheet1"){
		    my ($row_min, $row_max) = $worksheet->row_range();
		    my ($col_min, $col_max) = $worksheet->col_range();
		    
		    for my $row ($row_min..$row_max){
			my @row_values;
			for my $col ($col_min..$col_max){
			    
			    my $cell = $worksheet->get_cell($row, $col);
			    
			    if($cell){
				push(@row_values, $cell->value());
			    }
			}
			
			my $row_entry = join("\t", @row_values);
			push(@excel_output, $row_entry);
			    
			
		    }
		}
	    }
	    
	    open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
	    foreach my $entry (@excel_output){
		
		print OUTFILE $entry . "\n";
# 		die $entry . "\n";
	    }
	    close(OUTFILE) or die "Couldn't close file $outfile";
	}
    }
#    my $xls_files = find_files($input_dir, "xls");
#    foreach my $xls_file (sort keys %{$xls_files}){
#        my $infile = $xls_files->{$xls_file};
#        my $filename = $xls_file;
#        $filename =~ s/vs/v/g;
#        $filename =~ s/v/vs/g;
#        $filename =~ s/\.xls/\.txt/g;
#        
#        my $outfile = join('/', $outfile_dir, $filename);
#        unless(-s $outfile){
#            print join(" ", "Parsing xls file:", $xls_file . "\n");
#            open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
#            open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
#            while(<INFILE>){
#                chomp $_;
#                $_ =~ s/\r/\n/g;
#                
#                print OUTFILE $_ . "\n";
#            }
#            close(INFILE) or die "Couldn't close file $infile";
#            close(OUTFILE) or die "Couldn't close file $outfile";
#        }
#    }
    return $outfile_dir;
}

sub parse_txt_files{
    
    my $input_dir = shift;
	die "Error lost input directory" unless defined $input_dir;
    
	my $output_dir = shift;
	die "Error lost output directory" unless defined $output_dir;
    
    my $txt_files = find_files($input_dir, "txt");
    
    my (%sample_summary, %name_ids);
    my @sample_names;
    foreach my $txt_file (sort {$a cmp $b} keys %{$txt_files}){
	print join(" ", "Parsing txt file:", $txt_file . "\n");
	my $filename = $txt_files->{$txt_file};
	
	my $sample_name = $txt_file;
	$sample_name =~ s/\.txt//g;
	$sample_name=~ s/-/\./g;

	push(@sample_names, $sample_name);
	
	open(INFILE, "<$filename") or die "Couldn't open file $filename for reading, $!";
	my $i = 0;
	my (%column_index, %column_values);
	while(<INFILE>){
	    chomp $_;
	    $_ =~ s/ /\t/g;
	    my @split_entry = split(/\t/, $_);
	    if($i eq 0){
		my $j = 0;
		foreach my $column_name (@split_entry){
		    $column_index{$j} = $column_name;
		    $j++;
		}
	    }
	    else{
		my $j = 0;
		foreach my $column_value (@split_entry){
		    push(@{$column_values{$column_index{$j}}}, $column_value);
		    $j++;
		}
	    }
	    $i++;
	}
	close(INFILE) or die "Couldn't close file $filename";
	my $sizeFC = scalar @{$column_values{"FC"}};
	my $sizeAdjPVal = scalar @{$column_values{"adj.P.Val"}};
	warn "$sizeFC == $sizeAdjPVal\n";
	if($sizeFC eq $sizeAdjPVal){
	    my $size = ($sizeFC + $sizeAdjPVal)/2;
	    for(my $i = 0; $i < $size; $i++){
		my ($name, $id, $FC_value, $adjP_value) = (@{$column_values{"Name"}}[$i], @{$column_values{"ID"}}[$i], @{$column_values{"FC"}}[$i], @{$column_values{"adj.P.Val"}}[$i]);
		my $name_id = join("=", $name, $id);
# 		my $name_id = $id;
		$name_ids{$name_id} = $name_id;

		my $FC;
		if(($FC_value eq "N/A" or $FC_value eq "NA" or $FC_value eq "#VALUE!") 
		  and ($adjP_value eq "N/A" or $adjP_value eq "NA" or $adjP_value eq "#VALUE!")){
		  $FC = "null";
		  $sample_summary{$name_id}{$sample_name} = $FC;
		  next;
		}
		
		if($FC_value > 1.6 and $adjP_value < 0.05){
		    $FC = sprintf("%.9f", log_n($FC_value, 2));
		}
		elsif($FC_value < 0.6 and $adjP_value < 0.05){
		    $FC = sprintf("%.9f", log_n($FC_value, 2));
		}
		else{
		    $FC = "null";
		}
		$sample_summary{$name_id}{$sample_name} = $FC;
	    }
	}else{
	    die "sizes of FC and adj.P.Val are not equal";
	}
    }

    for my $sample_name (sort {$a cmp $b} @sample_names) {
	for my $name_id (sort {$a cmp $b} keys %name_ids) {
	
	  $sample_summary{$name_id}{$sample_name} = "null" unless(defined($sample_summary{$name_id}{$sample_name}));
	
	}
    }

    my %submission_summary;
    for my $name_id (sort keys %sample_summary){
	
	my @summary_row;
	for my $sample_name (sort keys %{$sample_summary{$name_id}}) {
#            push(@summary_row, "$sample_name ==> $sample_summary{$name_id}{$sample_name}");
	    push(@summary_row, $sample_summary{$name_id}{$sample_name});

	}
	$submission_summary{$name_id} = join("\t", $name_id, @summary_row);
    }
    
    my $outfile = join("/", $output_dir, "submission-summary.txt");
    open(OUTFILE1, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
    print OUTFILE1 join("\t", "ID_REF", @sample_names) . "\n";


  my $microarray_ids_file = join("/", $output_dir, "microarray-ids.txt");
  open(OUTFILE2, ">$microarray_ids_file") or die "Couldn't open file $microarray_ids_file for writting, $!";
    foreach my $name_id (sort {$a cmp $b} keys %submission_summary){
	
	    my @split_entry = split(/\t/, $submission_summary{$name_id});
	    my $is_null = "true";
	    for (my $j = 1; $j < scalar(@split_entry); $j++){
		  if($split_entry[$j] ne "null"){
			$is_null = "false";
			last;
		  }
	    }
	    if($is_null eq "false"){
		  print OUTFILE1 $submission_summary{$name_id} . "\n";
		  print OUTFILE2 $name_id . "\n";
	    }
	
    }
    close(OUTFILE1) or die "Couldn't close file $outfile";
    
#   close(OUTFILE2) or die "Couldn't close file $microarray_ids_file";
#   foreach my $name_id (sort {$a cmp $b} keys %name_ids){
#       
#       print OUTFILE $name_id . "\n";
#   }
#   close(OUTFILE) or die "Couldn't close file $microarray_ids_file";
}

sub log_n {
      my $number = shift;
      my $base = shift;
      return log($number)/log($base);
}