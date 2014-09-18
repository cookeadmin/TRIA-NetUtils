=head1 NAME
 
 B<extract_ipod_music.pl> - Extract music from your iPod.
 
 =head1 DESCRIPTION
 
 This program copies all the music from your iPod to your desktop
 or directory of choice using the MP4::Info Perl module found on CPAN.
 Download using the following command 'sudo ./install.sh -v' as the program
 is dependant on several modules. (-v is for VERBOSE)
 
=cut

#!/usr/bin/perl
use warnings;
use strict;

use Mail::Mailer;
use Switch;
use POSIX;

use constant DATE => strftime("%Y-%m-%d", localtime);
use constant TIME => strftime("%%H-%M-%S", localtime);

use Getopt::Long;
my ($input_dir, $output_dir);
my @options = (
    'i=s'    => \$input_dir,
    'o=s'    => \$output_dir,
);
&GetOptions(@options);

usage() unless (
    defined $input_dir
    and defined $output_dir
);

#my ($to_address,$from_address, $subject,$message) = ('kevin5@ualberta.ca','cookeadmin <root@cookelab.ccis.ualberta.ca>', "Testing Mail", "This is a test! The time is " . DATETIME . ". If you see this message the test passed");
my ($to_address,$from_address,$subject,$message) = ('14039195832@pcs.rogers.com','cookeadmin <root@cookelab.ccis.ualberta.ca>', "Testing Mail", join(" ", ""));

send_mail($to_address,$from_address,$subject,$message);

sub usage {
    
    die <<"USAGE";
    
Usage: $0 -i input_dir -o output_dir
    
    e.g. perl extract_ipod_music.pl -i /Volumes/YOURNAMEIPOD -o /Users/yourname/Desktop/MYIPOD_RECOVERED
    
    Description - This program copies all the music from your iPod to your desktop or directory of choice.
    
OPTIONS:
    -i input_dir - The volume drive that the iPod is
    located once plugged in using the usb drive on the
    computer.
    e.g. /Volumes/YOURNAMEIPOD
    
    -o output_dir - The output directory where you want the
    organized *.m4a containing artist folders to be
    stored.
    e.g. /Users/yourname/Desktop/MYIPOD_RECOVERED
    
Default: The desktop directory for your
    convinence.
    e.g. /Users/yourname/Desktop/YOURNAMEIPOD-MM-DD-YYYY
    
USAGE
}

# Make the output directory if it does not currently exist.
unless(-d $output_dir){
    mkdir $output_dir or die "Cannot create directory $output_dir: $!";
}


# Find all music files in the iPod music folder with the extension *.m4a.
# The input directory is /Volumes/YOURNAMEIPOD/iPod_Control/Music.
#my $ipod_music_dir = join('/', $input_dir, "iPod_Control", "Music");
#my ($music_files, $file_count) = find_music_files($ipod_music_dir);
#
## Contains all information on a file (original and new paths) and error status.
#my %copy_error_info;
#
## Iterate through the files with the extension *.m4a, rename them to the following
## format artist - title.m4a, and copy to the output folder $output_dir/$artist_dir/
#my $num_of_files = 0;
#foreach my $file (sort keys %{$music_files}){
#    # warn $music_files->{$file} , "\n";
#    
#    # Get the file has mp4 tag metadata.
#    my $tag = get_mp4tag($music_files->{$file});
#    
#    # Artist name directory and new music filename.
#    my($artist_dir, $music_file_name);
#    
#    # Checks if the file possesses any mp4tag metadata.
#    if(defined($tag)){ # Process mp4 tag metadata to generate new filename and artist directory.
#        my $mp4 = new MP4::Info $music_files->{$file};
#        my $artist = $mp4->artist;
#        my $title = $mp4->title;
#        $title =~ s/\/\s*//g; # Make sure '/' characters are filtered out.
#        $artist_dir = join('/', $output_dir, $artist);
#        $music_file_name = join(" - ", $artist, $title . ".m4a");
#    }
#    else{ # File does not possess any metadata to process.
#        $artist_dir = join('/', $output_dir, "untitled");
#        $music_file_name = $file;
#    }
#    
#    # Make the artist output directory if it does not currently exist.
#    unless(-d $artist_dir){
#        mkdir $artist_dir or die "Cannot create directory $artist_dir: $!";
#    }
#    
#    # Create new filename and copy to the specified location.
#    my $new_file_name = join("/", $artist_dir, $music_file_name);
#    warn join(" ", "Copying", $music_file_name, "to\n", $artist_dir) , "\n";
#    my $copy_status = copy("$music_files->{$file}","$new_file_name"); #or die "Copy failed: $!";
#    
#    # Checks the copy status of the file.
#    if($copy_status eq 1){ # Copy successful, count number of files copied successfully.
#        $num_of_files++;
#    }
#    else{ # Copy unsuccessful, obtain filename (original and new file paths) and error status.
#        $copy_error_info{$file}{"m4aEncoded"} = $music_files->{$file};
#        $copy_error_info{$file}{"m4aDecoded"} = $new_file_name;
#        $copy_error_info{$file}{"errMsg"} = "Copy failed: $!";
#    }
#}
#
#warn "$num_of_files out of $file_count .m4a music files copied successfully....\n";

## Tries to re-copy files that did not download properly.
#warn "File(s) with copy errors:\n";
#foreach my $file (sort keys %copy_error_info){
#    warn join("\n", $copy_error_info{$file}{"m4aEncoded"}, $copy_error_info{$file}{"errMsg"}, "Attempting to redo copying process....") . "\n";
#    my $copy_status = copy($copy_error_info{$file}{"m4aEncoded"}, $copy_error_info{$file}{"m4aDecoded"});
#    
#    if($copy_status eq 1){# If the copying process is successful, a message indicating that status appears.
#        warn join(" ", "File:", $copy_error_info{$file}{"m4aDecoded"}) , "\n", "copied successfully....\n";
#    }
#    else{ # If the copying process is unsuccessful, a message indicating that status appears, along with the
#        # original location of the file on the iPod and the name of the file so you can change it if
#        # you can manually copy the file required.
#        warn join(" ", "File:", $copy_error_info{$file}{"m4aDecoded"}) , "\n", join(" ", "Failed to copy from", $copy_error_info{$file}{"m4aEncoded"}, "successfully....") . "\n";
#    }
#}

=head1 SUBROUTINES
 
 (\%music_files, $file_count) = find_music_files($ipod_music_dir) - Find all music files in the specified iPod music folder with the extension *.m4a.
 
 Input paramater(s):
 
 $ipod_music_dir - iPod music input directory.
 
 Output paramater(s):
 
 \%music_files - A hash reference containing all the files with file extension *.m4a in key/value pairs.
 
 key => filename [A-Z]{4}.m4a ( e.g. ABCD.m4a )
 value => absolue filepath /path/to/[A-Z]{4} ( e.g. /path/to/ABCD.m4a )
 
 $file_count - The number of files stored with file extension *.m4a.
 
=cut
#sub find_structure_files{
#    
#    my $structure_dir = shift or die "lost structure dir";
#    
#    my $suffix = shift or die "lost file suffix";
#    
#    my (%structure_files, $file_count);
#    
#    $file_count = 0;
#    opendir(DIR, $structure_dir) || die "Error in opening dir $$structure_dir\n";
#    while( my $file_name  = readdir(DIR)){
#            my $structure_file_name = join('/', $structure_dir, $file_name);
#            warn "$structure_file_name\n" if ($file_name =~ m/\.$suffix$/);
#            $structure_files{$file_name} = $structure_file_name if($file_name =~ m/\.$suffix$/);
#            $file_count++ if($file_name =~ m/\.$suffix$/);
#        }
#    }
#    closedir(DIR);
#    return (\%structure_files, $file_count);
#}

=pod
 
 $output_dir = get_output_dir($input_dir) - Generate iPod output directory using the specified iPod input directory.
 
 Input paramater(s):
 
 $input_dir - The iPod input directory.
 
 Output paramater(s):
 
 $output_dir - The iPod output directory.
 
=cut
sub send_mail{
    my $to_address = shift or die "lost Recipient email address";
    my $from_address = shift or die "lost Sender email address";
    my $subject = shift or die "lost email subject";
    my $message = shift or die "lost email message";
    
    warn "Sending mail with $subject to $to_address\n";
    my $mailer = new Mail::Mailer("sendmail");
    $mailer->open( {
        To       => $to_address,
        From     => $from_address,
        Subject  => $subject
    } );
    
    print $mailer <<END_OF_MESSAGE;
    $message
END_OF_MESSAGE
        
    close $mailer;
}

=head1 SYNOPSIS
 
 B<Usage:> extract_ipod_music.pl -i input_dir -o output_dir
 
 e.g. sudo perl extract_ipod_music.pl -i /Volumes/YOURNAMEIPOD -o /Users/yourname/MYIPOD_RECOVERED
 
 B<OPTIONS:>
 
 -i input_dir - The volume drive that the iPod is
 located once plugged in using the usb drive on the
 computer.
 e.g. /Volumes/YOURNAMEIPOD
 
 -o output_dir - The output directory where you want the
 organized *.m4a containing artist folders to be
 stored.
 e.g. /Users/yourname/Desktop/MYIPOD_RECOVERED
 
 B<Default:> The desktop directory for your
 convinence.
 e.g. /Users/yourname/Desktop/YOURNAMEIPOD-MM-DD-YYYY
 
 =head1 AUTHOR
 
 B<Kevin Muirhead 2014>
 
 GPL License (Open Source)
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but B<WITHOUT ANY WARRANTY>; without even the implied warranty of
 B<MERCHANTABILITY> or B<FITNESS FOR A PARTICULAR PURPOSE>.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see L<http://www.gnu.org/licenses/>.
 
=cut