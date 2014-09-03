#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;

# perl xml_parser.pl -i ~/workspace/adriana/kegg_pathway/ko00902.xml -o ~/workspace/adriana/kegg_pathway

my ($infile, $output_dir);
GetOptions(
      'i=s'    => \$infile,
      'o=s'    => \$output_dir,
);
usage() unless (
      defined $infile
      and defined $output_dir
);


sub usage {
    
die <<"USAGE";
    
Usage: $0 -i infile -o output_dir
    
Description - 
    
OPTIONS:
      -i infile - 

      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}




my $infilename = fileparse($infile);

my @xml_content = ();
open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
while(<INFILE>){
      chomp $_;
#       warn "$_\n";
      push(@xml_content, $_);
}
close(INFILE) or die "Couldn't close file $infile";

my $xml_data = join("\n", @xml_content) . "\n";
$xml_data =~ s/\n\s+type=/ type=/g;

$xml_data =~ s/"//g;

my @xml_formated_data = split(/\n/, $xml_data);

my $filename = fileparse($infile, qr/\.[^.]*/);
my $outfile = join('/', $output_dir, $filename . ".tsv");
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE join("\t", "name","x_coord","y_coord","width","height") . "\n";
foreach my $data (@xml_formated_data){

      if($data =~ m/<graphics name=(.+) fgcolor=#\w+ bgcolor=#\w+ type=rectangle x=(\d+) y=(\d+) width=(\d+) height=(\d+)\/>/){
	    my ($name,$x_coord,$y_coord,$width,$height) = ($1,$2,$3,$4,$5);
	    print OUTFILE join("\t", $name,$x_coord,$y_coord,$width,$height) . "\n";
      }

}
close(OUTFILE) or die "Couldn't close file $outfile"; 