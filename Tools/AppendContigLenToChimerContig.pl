#!/bin/perl

# Append the contig lengths to the chimerContig file
# usage: a.pl ref.fa xxx.chimerContig >...

use strict ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my %length ;
my $seq = "" ;
my $id ;
while (<FP1>)
{
	chomp ;
	if ( /^>/ )
	{
		if ( $seq ne "" )
		{
			$length{$id} = length( $seq ) ;
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
	$length{$id} = length( $seq ) ;
}

my @cols ;
while ( <FP2> )
{
	chomp ;
	@cols = split /\s+/ ;
	#print $cols[0], "\n" ;
	print $_, "\t", $length{">".$cols[0]}, "\t", $length{">".$cols[2]}, "\n" ;
}