#!/bin/perl

# Reformat the coords file to the format of ChainNucmer
# Usage: a.pl origin.oneline.fa xxx.coords [-reverse|-self]

use strict ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my $reverse = 1 ;
my $self = 0 ;
if ( $ARGV[2] eq "-reverse" ) # Notice the definition of this parameter is not straightforward
{
	$reverse = 0 ;
}
elsif ( $ARGV[2] eq "-self" )
{
	$self = 1 ;
}

my %origLength ;
my %output ;
my %fragCnt ;
my @cols ;
my $text ;
my $line ;
my $i ;
while ( <FP1> ) # Assume one-line fasta format
{
	chomp ;
	my $id = (split /\s+/, substr( $_, 1 ))[0] ; 
	my $seq = <FP1> ;
	chomp $seq ;
	@cols = split /N+/, $seq ;
	my $len = 0 ;
	for ( $i = 0 ; $i < scalar( @cols ) ; ++$i )
	{
		$len += length( $cols[$i] ) ;
	}
	# The length without gap
	#$origLength{$id} = length( $seq ) ;
	$origLength{$id} = $len ;
	$fragCnt{$id} = 0 ;
}

#     3939     4038  |       25      124  |      100      100  |    99.00  | scaffold_202.56      scf620008933730
$line = <FP2> ;
$line = <FP2> ;
$line = <FP2> ;
$line = <FP2> ;
$line = <FP2> ;
while ( <FP2> )
{
	chomp ;
	$line = $_ ;
	@cols = split /\s+/, $line ;
	shift @cols if ( $cols[0] eq "" ) ;

	my $start ; 
	my $end  ;
	my $destId ; 
	my $dstart ;
	my $dend  ;
	my $id ;

	if ( $reverse == 0 )
	{
		$start = $cols[0] ;
		$end = $cols[1] ;
		$destId = $cols[12] ;
		$dstart = $cols[3] ;
		$dend = $cols[4] ;
		$id = $cols[11] ;
	}
	else
	{
		$start = $cols[3] ;
		$end = $cols[4] ;
		$destId = $cols[11] ;
		$dstart = $cols[0] ;
		$dend = $cols[1] ;
		$id = $cols[12] ;

		if ( $start > $end )
		{
			$start = $cols[4] ;
			$end = $cols[3] ;
			$dstart = $cols[1] ;
			$dend = $cols[0] ;
		}
	}
	$text = "$start $end $destId $dstart $dend\n" ;
	$output{ $id } .= $text ;
	++$fragCnt{$id} ;

	# Add the other way around
	if ( $self == 1 )
	{
		($start, $end, $dstart, $dend) = ($dstart, $dend, $start, $end) ;
		if ( $start > $end )
		{
			($start, $end) = ( $end, $start ) ;
			($dstart, $dend) = ($dend, $dstart ) ;
		}
		($id, $destId) = ($destId, $id) ;
		$text = "$start $end $destId $dstart $dend\n" ;
		$output{ $id } .= $text ;
		++$fragCnt{$id} ;
	}
}

foreach my $key (keys %origLength )
{
	print "$key $origLength{$key} $fragCnt{$key}\n$output{$key}" ;
}
