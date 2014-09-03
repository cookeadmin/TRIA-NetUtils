
#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use File::Basename;
use WWW::Mechanize;
use open qw(:std :utf8);

# perl extract_html_tables.pl -i /home/cookeadmin/workspace/adriana/kegg_pathway -o /home/cookeadmin/workspace/adriana/kegg_pathway

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
my $mech = WWW::Mechanize->new;
open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
my $i = 0;
my %html_links = ();
while(<INFILE>){
      chomp $_;
      if($i ne 0){
	    my @split_entry = split(/\t/, $_);
	    my ($name,$html_source) = @split_entry;
	    warn $_ . "\n";

	    push(@{$html_links{$name}}, $html_source);
# 	    my $outfile = join("/", $output_dir, $name . ".png");
# 	    system("convert /home/cookeadmin/workspace/adriana/kegg_pathway/ko00902.png -stroke black -fill \"rgba( 255, 215, 0 , 1.0 )\" -draw \"rectangle $x_coord,$y_coord,$width,$hieght \" $outfile");
      }
      $i++;
}
close(INFILE) or die "Couldn't close file $infile";

foreach my $gene_name (sort keys %html_links){
      warn "Processing gene id: " . $gene_name . ".....\n";
      my $outfile = join("/", $output_dir, $gene_name . "-genes.html");
      open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
      foreach my $html_source (@{$html_links{$gene_name}}){
	    
	    warn $html_source . "\n";
	    my $response = $mech->get($html_source);
	    print OUTFILE $mech->content;
	    sleep(5); 
      }
      close(OUTFILE) or die "Couldn't close file $outfile";
}
     