#!/bin/perl

# Verify regions in the input file from other evidences
# usage: a.pl fff.oneline.fasta xxx.region [options] > output
# The first four column of the region file is:
# scaf_id start end length ...

# options:
# 	-d a.depth
#   -iv b.intervalVerify

use strict ;

my @regions ;
my %scafLength ;
my @cols ;
my %falseRegion ;

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
	$scafLength{$id} = length( $seq ) ;
}

while ( <FP2> )
{
	chomp ;
	@cols = split ;
	my $text = $cols[0]." ".$cols[1]." ".$cols[2] ;
	push @regions, $text ;
}

sub VerifyDepth
{
	open FP3, $_[0] ;
	
	my %depth ;
	my $id ;
	my $i ;
	my $j ;
	my @cols ;
	
	for $id (keys %scafLength )
	{
		@#{$depth{$id}} = $scafLength{$id} + 1 ;
	}
	
	while ( <FP3> )
	{
		chomp ;
		@cols = split ;
		@{$depth{$cols[0]}}[$cols[1]] = $cols[2] ;
	}
	
	for	( $i = 0 ; $i < scalar( @regions ) ; ++$i )
	{
		@cols = split /\s+/, $regions[$i] ;
		my $array = \@{$depth{$cols[0]}} ;
		for ( $j = $cols[1] ; $j < $cols[2] ; ++$j )
		{
			if ( $array->[$j] < 1 ) 
			{
				$falseRegion{ $regions[$i] } .= "Hallow from $_[0]. " ;
				last ;
			}
		}
	}
	close FP3 ;
}

sub VerifyIntervalVerify
{
	open FP3, $_[0] ;
	
	my %verification ;
	my $id ;
	my $i ;
	my $j ;
	my @cols ;
	my $line ;
	
	while ( <FP3> )
	{
		chomp ;
		$id = (split /\s+/, substr( $_, 1 ))[0] ; 
		$line = <FP3> ;
		chomp $line ;
		$verification{ $id } = $id ;
	}
	
	for	( $i = 0 ; $i < scalar( @regions ) ; ++$i )
	{
		@cols = split /\s+/, $regions[$i] ;
		if ( substr( $verification{ $cols[0] }, $cols[1], $cols[2] - $cols[1] + 1 ) =~ /0/ ) 
		{
			$falseRegion{ $regions[$i] } .= "Unconfirmted portion from $_[0]. " ;
			next ;
		}
	}	
}


my $argvInd ;
for ( $argvInd = 2 ; $argvInd < scalar( @ARGV ) ; ++$argvInd )
{
	if ( $ARGV[ $argvInd ] eq "-d" )
	{
		VerifyDepth( $ARGV[$argvInd + 1]) ;
		++$argvInd ;
	}
	elsif ( $ARGV[ $argvInd ] eq "-iv" )
	{
		VerifyIntervalVerify( $ARGV[ $argvInd + 1] )  ;
		++$argvInd ;
	}
}

# Output the regions that can not be verified
my $key ;
for $key (keys %falseRegion)
{
	@cols = split /\s+/, $key ;
	print "$key ".($cols[2] - $cols[1] + 1)." $falseRegion{$key}\n" ;
}
 