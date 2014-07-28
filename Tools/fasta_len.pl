#!/bin/perl

my @scafLength ;
my $seq = "" ;

while (<>)
{
	if ( /^>/ )
	{
		if ( $seq ne "" )
		{
			push @scafLength, length( $seq ) ;
		}
		$seq = "" ;
	}
	else
	{
		chomp ;
		$seq .= $_ ;
	}
}
push @scafLength, length( $seq ) ;

my $i = 0 ;
my $culm = 0 ;
foreach my $l (sort {$b<=>$a} @scafLength)
{
	$culm += $l ;
	print $i, "\t", $l, "\t", $culm, "\n" ;
	++$i ;
}