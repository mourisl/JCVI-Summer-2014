#!/bin/perl

# Convert scaffolds file to contig file
# Usage: a.pl xxx.fa N_size  > yy.fa
use strict ;

my $n ;
$n = pop @ARGV ;
my $seq = "" ;
my $id = "" ;
my @cols ;
my $i ;

while (<>)
{
	chomp ;
	if ( /^>/ )
	{
		if ( $seq ne "" )
		{
			@cols = split /[N]{$n,}/, $seq ;
			for ( $i = 0 ; $i < scalar( @cols ) ; ++$i )
			{
				print $id, "_$i\n", $cols[$i], "\n" ;
			}
		}
		$id = $_ ;
		$seq = "" ;
	}
	else
	{
		$seq .= $_ ;
	}
}
if ( $seq ne "" )
{
	@cols = split /[N]{$n,}/, $seq ;
	for ( $i = 0 ; $i < scalar( @cols ) ; ++$i )
	{
		print $id, "_$i\n", $cols[$i], "\n" ;
	}
}	