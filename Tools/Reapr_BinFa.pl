#!/bin/perl

# Compare the bin.fa with the .fa file to see the modifications within a contig.
# usage: a.pl bin.fa .fa

use strict ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my %scaf ;
my $seq ;
my $scafId = "" ;
my @cols ;
my @cols2 ;

while (<FP2>)
{
	if ( /^>/ )
	{
		if ( $scafId ne "" )
		{
			$scaf{ $scafId } = $seq ;
		}
		
		$seq = "" ;
		chomp ;
		$scafId = substr( $_, 1 ) ;
	}
	else
	{
		chomp ;
		$seq = $seq.$_ ;
	}
}
if ( $scafId ne "" )
{
	$scaf{ $scafId } = $seq ;
}

# Process the bin.fa file
$scafId = "" ;
while ( <FP1> )
{
	if ( /^>/ )
	{
		if ( $scafId ne "" )
		{
			my $key ;	
			my $val ;
			@cols2 = split /_/, $scafId ;
			#print $cols2[1], "\n" ;
			while ( ($key, $val) = each( %scaf ) )
			{
				@cols = split /_/, $key ;
				#print $cols[1], "\n" ;
				next if ( $cols[0] ne $cols2[0] || $cols[1] > $cols2[1] || $cols[2] < $cols2[2] ) ;
				print $scafId, "\n", substr( $scaf{ $key }, $cols2[1] - $cols[1], $cols2[2] - $cols2[1] + 1 ), "\n" ;
				#last ;
			}
		}
		$seq = "" ;
		chomp ;
		$scafId = substr( $_, 11 ) ;
	}
	else
	{
		chomp ;
		$seq = $seq.$_ ;
	}
}
if ( $scafId ne "" )
{
	my $key ;	
	my $val ;
	@cols2 = split /_/, $scafId ;
	while ( ($key, $val) = each( %scaf ) )
	{
		@cols = split /_/, $key ;
		
		next if ( $cols[0] ne $cols2[0] || $cols[1] > $cols2[1] || $cols[2] < $cols2[2] ) ;
		print $scafId, "\n", substr( $scaf{ $key }, $cols2[1] - $cols[1], $cols2[2] - $cols2[1] + 1 ), "\n" ;			
		#last ;
	}
}
