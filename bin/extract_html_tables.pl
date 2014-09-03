#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use HTML::TableExtract;
use open qw(:std :utf8);

# perl extract_html_tables.pl -i /home/cookeadmin/workspace/adriana/kegg_pathway -o /home/cookeadmin/workspace/adriana/kegg_pathway

my ($input_dir, $output_dir);
GetOptions(
      'i=s'    => \$input_dir,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $input_dir
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -i input_dir -o output_dir
    
Description - 
    
OPTIONS:
      -i input_dir - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my $files = find_files($input_dir, "html");

foreach my $file (keys %{$files}){
      my $infile = $files->{$file};

      my @html_content = ();
      open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
      while(<INFILE>){
	    chomp $_;
      #       warn "$_\n";
	    push(@html_content, $_);
      }
      close(INFILE) or die "Couldn't close file $infile";

      my $html_data = join("\n", @html_content) . "\n";


      my $te = HTML::TableExtract->new();
      $te->parse($html_data);

      my $filename = fileparse($infile, qr/\.[^.]*/);
      my $outfile = join('/', $output_dir, $filename . ".txt");
      open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
      foreach my $ts ($te->tables()){
	    foreach my $row ($ts->rows()){
		  if($row){
			my $entry = join("\t", @$row );
			$entry =~ s/\t//g;
			print OUTFILE $entry . "\n";
		  }
	    }
      }
      close(OUTFILE) or die "Couldn't close file $outfile"; 
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