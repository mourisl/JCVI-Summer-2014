#!/bin/perl

# Evaluate the verified bases and the new N50 from the IntervalVerify file

# usage: a.pl [-breakN aaa.fa N_size] < xxx.IntervalVerify

use strict ;

my $i ;
my @newScaf ;
my $verified ;
my $sum = 0 ;
my $oldSum = 0 ;
my @cols ;
my $id ;
my $N50Sum = 0 ;
my $breakN = 0;
my $breakNSize = 20 ;
my %scafSeq ;

for ( $i = 0 ; $i < scalar( @ARGV ) ; ++$i )
{
	if ( $ARGV[$i] eq "-breakN" )
	{
		$breakN = 1 ;
		$breakNSize = $ARGV[$i + 2] ;
		open FP1, $ARGV[$i + 1] ;
		while ( <FP1> )
		{
			chomp ;
			if ( /^>/ )
			{
				$id = $_ ;
			}
			else
			{
				$scafSeq{ $id } .= $_ ;
			}
		}
		$i += 2 ;
	}
}
undef @ARGV ;

while ( <STDIN> )
{
	chomp ;
	if ( /^>/ )
	{
		$id = $_ ;
		next ;
	}
	$verified = $_ ;
	
	if ( $breakN eq 1 )
	{
		# Set to 0 if it is in a stretch of N
		while ( $scafSeq{$id} =~ /N+/g )
		{
			if ( $+[0] - $-[0] >= $breakNSize )
			{
				my $j ;
				for ( $j = $+[0] ; $j < $+[1] ; ++$j )
				{
					substr( $verified, $j, 1, '0' ) ;
				}
			}
		}
	}
	
	@cols = split /0/, $verified ;
	$oldSum += length( $verified ) ;
	for ( $i = 0 ; $i < scalar( @cols ) ; ++$i )
	{
		push @newScaf, length( $cols[ $i ] ) ; 
	}

	if ( scalar( @cols ) == 0 )
	{
		print "Dubious scaffold: $id\n" ;
	}
}

for ( $i = 0 ; $i < scalar( @newScaf ) ; ++$i )
{
	$sum += $newScaf[$i] ;
	$N50Sum += $newScaf[$i] if ( $newScaf[$i] >= 2000 ) ;
}
print "Total bases: $oldSum\nVerified bases: $sum\nVerified bases from stretch no less than 2000bp: $N50Sum\n" ;

my $psum = 0 ;
my $k = 10 ;
#$N50Sum = 192921369 ;

foreach $i (sort {$b<=>$a} @newScaf )
{
	last if ( $i < 2000 ) ;
	$psum += $i ;
	while ( $psum >= $N50Sum * ( $k / 100 ) )
	{
		print "N$k = $i\n" ;
		$k += 10 ;
	}
}
