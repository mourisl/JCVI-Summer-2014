#!/bin/perl

# Build the histo files for each type of coverage.
# Usage: a.pl bin_size [-cumulative]< xxx.chain > out
# scaffold_908.0 122 122 100.000000 scaffold0572 1 [(1 122) (10238 10119)]

use strict ;

my @coveredLength ;
my @totalLength ;
my @number ;

my $bin = $ARGV[0] ;
my $i ;
my @cols ;

for ( $i = 0 ; $i <= 100 ; $i += $bin )
{
	$coveredLength[ $i / $bin ] = 0 ;
	$totalLength[ $i / $bin ] = 0 ;
	$number[ $i / $bin ] = 0 ;
}

while ( <stdin> )
{
	chomp ;
	@cols = split ;
	my $ind = $cols[3] / $bin ;

	$coveredLength[ $ind ] += $cols[1] ;
	$totalLength[ $ind ] += $cols[2] ;
	++$number[$ind] ;
}

for ( $i = 0 ; $i <= 100 ; $i += $bin )
{
	print "$i $coveredLength[$i/$bin] $totalLength[$i / $bin] $number[$i / $bin]\n" ;
}
