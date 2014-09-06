#!/usr/bin/perl
use warnings;
use strict;

my $num_chars = 10;
my $pwd;
my @my_char_list = (('A'..'Z'), ('a'..'z'), ('!','@','%','^'), (0..9));
my $range_dis = $#my_char_list + 1;
for (1..$num_chars) {
   $pwd .= $my_char_list[int(rand($range_dis))];
}
print "$pwd\n";