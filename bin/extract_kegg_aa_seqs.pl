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

my $files = find_files($input_dir, "txt");
my $plantdb_codes = get_plantdb_codes();
my %html_links = ();
foreach my $file (keys %{$files}){
      warn $file . "\n";
      my $infile = $files->{$file};
      my $gene_name = fileparse($infile, qr/\.[^.]*/);
      $gene_name =~ s/-genes//g;
      my @html_content = ();
      open(INFILE, "<$infile") or die "Couldn't open file $infile for reading, $!";
      while(<INFILE>){
	    chomp $_;
# 	    warn $_ . "\n";
	    push(@html_content, $_);
      }
      close(INFILE) or die "Couldn't close file $infile";
      my @split_html_content = split(/Download RDF   /, join("\n", @html_content) . "\n");
      my @aa_sequences = ();
      foreach my $aa_entry (@split_html_content){
	    $aa_entry =~ s/\n\n/\n/g;
	    $aa_entry =~ s/\n{3,}/\n/g;
# 	    warn $aa_entry . "\n";
	    my ($gene_db_name, $gene_kegg_name, $org_code, $org_name, $aa_length);
	    my @aa_seq = ();
	    my $is_aa_seq_found = "false";
	    my @split_aa_entry = split(/\n/, $aa_entry);
	    foreach my $entry (@split_aa_entry){
		  if($entry =~ /(.+)CDS.+/){
			$gene_db_name = $1;
			$gene_db_name =~ s/\s+//g;
# 			warn $gene_db_name . "\n";
		  }
		  if($entry =~ /^(K\d+.+)/){
			$gene_kegg_name = $1;
			$gene_kegg_name =~ s/\s\s/ /g;
# 			warn $gene_kegg_name . "\n";
		  }
		  if($entry =~ /Organism(.+)/ and $entry !~ /Organismal Systems/){
			my $organism = $1;
			($org_code, $org_name) = split(/\s\s/, $organism);
			$org_name =~ s/ \(.+\)//g;
#  			warn $organism . "\n";
		  }

		  if($entry =~ /AA seq(\d+) aa/){
			$aa_length = $1;
			$is_aa_seq_found = "true";
		  }elsif($entry =~ m/([A-Z]+)/ and $entry !~ /NT seq\d+ nt/ and $is_aa_seq_found eq "true"){
			my $aa_partial = $1;
			push(@aa_seq, $aa_partial);
		  }
		  if($entry =~ /NT seq\d+ nt/ and $is_aa_seq_found eq "true"){
			
			$is_aa_seq_found = "false";
		  }

	    }
	    if(defined($gene_db_name) and defined($gene_kegg_name) and defined($org_code) and defined($org_name) and defined($aa_length)){
		  my ($aa_header, $aa_sequence, $fasta_entry) = "";
		  $aa_header = join(" ", ">$gene_kegg_name", join(":", $org_code,$gene_db_name), "[$org_name]", join(' ', "length=$aa_length", "aa"));
		  $aa_sequence = join("\n", @aa_seq);
		  $fasta_entry = join("\n", $aa_header, $aa_sequence);
		  push(@aa_sequences, $fasta_entry);
	    }else{
		  next;
	    }
      }
      my $outfile = join('/', $output_dir, join("-", $gene_name , "proteins.fasta"));
      open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
      foreach my $fasta_entry (@aa_sequences){
	    print OUTFILE $fasta_entry . "\n";
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