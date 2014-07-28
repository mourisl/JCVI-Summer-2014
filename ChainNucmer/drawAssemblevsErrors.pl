#!/bin/perl

#usage: a.pl < assemble.chain
# scaffold_596.2 12382 12382 100.000000 chr3 1 [(1 12382) (4571774 4584131)]

my @cols ;
my %verified ;
my %assembleLen ;
my @idArray ;
my $len ;
my $sum = 0 ;
my $sumError = 0 ;
my $i = 0 ;
while ( <stdin> )
{
	chomp ;
	@cols = split ;
	push @idArray, "$i" ;
	$assembleLen{ $i } = $cols[2] ;
	$verified{$i} = $cols[1] ;
	++$i ;
}

foreach $i (sort {$assembleLen{$b}<=>$assembleLen{$a}} (@idArray))
{
	$len = $assembleLen{ $i };
	$sum += $len ;
	#print "== $len $verified{$i}\n" ;
	$sumError += ( $len - $verified{ $i } );
	print "$sum $sumError\n" ;
}
