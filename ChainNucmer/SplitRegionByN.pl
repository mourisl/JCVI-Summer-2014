#!/bin/perl

#usage: a.pl xxx.oneline.fa yyy.region

my %scaf ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

while ( <FP1> ) # Assume one-line fasta format
{
	chomp ;
	my $id = (split /\s+/, substr( $_, 1 ))[0] ; 
	my $seq = <FP1> ;
	chomp $seq ;
	# The length without gap
	#$origLength{$id} = length( $seq ) ;
	$scaf{$id} = $seq ;
}

my @cols ;
while ( <FP2> )
{
	chomp ;
	@cols = split /\s+/, $_ ;
	my $text = substr( $scaf{$cols[0]}, $cols[1] - 1, $cols[3] ) ;
	
	my $i ;
	my $j ;
	my $start = -1 ;
	for ( $i = 0 ; $i < length( $text ) ; ++$i )
	{
		if ( substr( $text, $i, 1 ) eq "N" )
		{
			if ( $start != -1 )
			{
				print "$cols[0] ", $cols[1] + $start, " ", $cols[1] + $i - 1, " ", $i - $start ;
				for ( $j = 4 ; $j < scalar( @cols ) ; ++$j )
				{
					print " $cols[$j]" ;
				}
				print "\n" ;
			}
			$start = -1 ;
		}
		else
		{
			$start = $i if ( $start == -1 ) ;			
		}
	}
	if ( $start != -1 )
	{
		print "$cols[0] ", $cols[1] + $start, " ", $cols[1] + $i - 1, " ", $i - $start ;
 
		for ( $j = 4 ; $j < scalar( @cols ) ; ++$j )
		{
			print " $cols[$j]" ;
		}
		print "\n" ;	
	}
}