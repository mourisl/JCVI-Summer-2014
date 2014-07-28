#!/bin/perl

# Compute how many fragments covered a position in the assembly
# Usage: a.pl origin.fa xxx.coords similiarity [-reverse]

use strict ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my %origLength ;
my %output ;
my %fragCnt ;
my @cols ;
my $text ;
my $line ;
my $i ;
my %coverage ;
my %scaf ;
my $reverse = 0 ;

$reverse = 1 if ( $ARGV[3] eq "-reverse" ) ;

while ( <FP1> ) # Assume one-line fasta format
{
	chomp ;
	my $id = (split /\s+/, substr( $_, 1 ))[0] ; 
	my $seq = <FP1> ;
	chomp $seq ;
	# The length without gap
	#$origLength{$id} = length( $seq ) ;
	$origLength{$id} = length( $seq ) ;
	$scaf{$id} = $seq ;
	@#{ $coverage{$id} } = length( $seq ) + 1 ;
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
	#print "$cols[10]\n" ;
	next if ( $cols[9] < $ARGV[2] ) ;
	if ( $reverse == 0 )
	{
		my $array = \@{$coverage{ $cols[11] } };
		
		for ( $i = $cols[0] ; $i <= $cols[1] ; ++$i )
		{
			++( $array->[$i] ) ;
		}
	}
	else
	{
		my $array = \@{$coverage{ $cols[12] } };
		my $start ;
		my $end ;
		
		if ( $cols[3] < $cols[4] )
		{
			$start = $cols[3] ;
			$end = $cols[4] ;
		}
		else
		{
			$start = $cols[4] ;
			$end = $cols[3] ;			
		}
		
		for ( $i = $start ; $i <= $end ; ++$i )
		{
			++( $array->[$i] ) ;
		}	
	}
}

for my $id (keys %origLength)
{
	my $array = \@{$coverage{ $id } } ;
	my $seq = \$scaf{$id} ;
	my $output = 0 ;
	for ( $i = 1 ; $i <= $origLength{$id} ; ++$i )
	{
		if ( substr( $$seq, $i - 1, 1 ) eq 'N' )
		{
			print "$id $i N\n" ;
			$output = 1 ;
			next ;
		}
		if ( ( defined $array->[$i] ) && $array->[$i] > 0 )
		{
			print "$id $i $array->[$i]\n" ;
			$output = 1 ;
		}
	}
	if ( $output == 0 )
	{
		print "$id 0 E\n" ;
	}
}