#!/usr/bin/perl
use strict;
use warnings;

#To run:
#perl batch_raxml.pl listfile.txt
#where listfile is a text file listing the names of all multiple sequences alignment files in phylip format (one file name per line).

my $listFile = $ARGV[0];
my @list;

open (AFILE, $listFile) or die "cannot open $listFile\n";
while (my $line = <AFILE>) {
        chomp $line;
        push @list, $line;
}
close AFILE;
print "\n@list\n\n"; #to test the elements of the array

for (my $i=0; $i<@list; $i++) {
        my $file = $list[$i];
        system("raxmlHPC -s $file -n $file -m GTRGAMMA -x 12345 -# 100 -f a -T 8 -p 12345");
}
