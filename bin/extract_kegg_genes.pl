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


my $html_prefix = 'http://www.genome.jp/dbget-bin/www_bget?';

my $files = find_files($input_dir, "txt");
my $plantdb_codes = get_plantdb_codes();
my %html_links = ();
foreach my $file (keys %{$files}){
      warn $file . "\n";
      my $infile = $files->{$file};
      my $gene_name = fileparse($infile, qr/\.[^.]*/);

      open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
      while(<INFILE>){
	    chomp $_;
	    if($_ =~ m/^([A-Z]+):(.+)/){
# 		warn "$_\n";
		  my ($dbname, $gene_names) = ($1, $2);
		  $dbname = lc($dbname);
		  $gene_names =~ s/\s//;
# 		  warn "$gene_name\n";
		  next unless(defined($plantdb_codes->{$dbname}));
		  my @split_gene_name_entries = split(/\s/, $gene_names);
		  my @gene_names = ();
		  foreach my $gene_name_entry (@split_gene_name_entries){
			$gene_name_entry =~ s/\(.+\)//;
		  	warn "$gene_name_entry\n";
			push(@gene_names, $gene_name_entry);
		  }
		  
		  my $gene_queries =  join("+", @gene_names);
		  my $html_url = join('', $html_prefix,$dbname,":",$gene_queries);
		  push(@{$html_links{$gene_name}}, $html_url);
	    }
      }
      close(INFILE) or die "Couldn't close file $infile";
}


my $outfile = join('/', $output_dir, "kegg_seqinfo_links.txt");
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE join("\t", "name","html_source") . "\n";
foreach my $gene_name (sort keys %html_links){
      foreach my $url (@{$html_links{$gene_name}}){
	    print OUTFILE join("\t", $gene_name, $url) . "\n";
      }
}
close(OUTFILE) or die "Couldn't close file $outfile"; 

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

sub get_plantdb_codes{
      my %plantdb_codes  = (
	    'ath' => 'Arabidopsis thaliana (thale cress)',
	    'aly' => 'Arabidopsis lyrata (lyrate rockcress) ',
	    'crb' => 'Capsella rubella',
	    'eus' => 'Eutrema salsugineum',
	    'cit' => 'Citrus sinensis (Valencia orange)',
	    'cic' => 'Citrus clementina (mandarin orange)',
	    'tcc' => 'Theobroma cacao (cacao)',
	    'gmx' => 'Glycine max (soybean)',
	    'pvu' => 'Phaseolus vulgaris (common bean)',
	    'mtr' => 'Medicago truncatula (barrel medic)',
	    'cam' => 'Cicer arietinum (chickpea)',
	    'fve' => 'Fragaria vesca (woodland strawberry)',
	    'pper' => 'Prunus persica (peach)',
	    'csv' => 'Cucumis sativus (cucumber)',
	    'rcu' => 'Ricinus communis (castor bean)',
	    'pop' => 'Populus trichocarpa (black cottonwood)',
	    'vvi' => 'Vitis vinifera (wine grape)',
	    'sly' => 'Solanum lycopersicum (tomato)',
	    'sot' => 'Solanum tuberosum (potato)',
	    'osa' => 'Oryza sativa japonica (Japanese rice) (RefSeq)',
	    'dosa' => 'Oryza sativa japonica (Japanese rice) (RAPDB)',
	    'obr' => 'Oryza brachyantha (malo sina)',
	    'bdi' => 'Brachypodium distachyon',
	    'sbi' => 'Sorghum bicolor (sorghum)',
	    'zma' => 'Zea mays (maize)',
	    'sita' => 'Setaria italica (foxtail millet)',
      );

      return \%plantdb_codes;
}