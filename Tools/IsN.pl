#!/bin/perl

# Read in a scaffold, output each position and whether it is N.
# Usage: perl a.pl < xxx.fa > yyy.out

my $line = <STDIN> ;
$line = <STDIN> ; # Suppose a one-line fasta.
my $i ;

for ( $i = 0 ; $i < length( $line ) ; ++$i )
{
	if ( substr( $line, $i, 1 ) eq "N" )
	{
		print "$i\t1\n" ; 
	}
	else
	{
		print "$i\t0\n" ;
	}
}