#!/bin/perl

# usage: a.pl < xxx.chain
# scaffold_596.2 12382 12382 100.000000 chr3 1 [(1 12382) (4571774 4584131)]
my @cols ;
my @aligned ;
my $i ;
my $len ;
my $sum ;

undef @aligned ;

while (<stdin>)
{
	#print $_ ;
	chomp ;
	@cols = split ;
	my $start ;
	my $end ;
	for ( $i = 6 ; $i < scalar( @cols ) ; $i += 4 ) 
	{
		if ( $cols[$i] =~ /([0-9]+)/ )
		{
			$start = $1 ;
		}
		else
		{
			die ;
		}
		
		if ( $cols[$i + 1] =~ /^([0-9]+)/ )
		{
			$end = $1 ;
		}
		else
		{
			die ;
		}		
		#print "$cols[$i] $cols[$i + 1] $start $end\n" ;
		push @aligned, $end - $start + 1 ;
	}
}

$i = 1 ;
$sum = 0 ;
foreach $len (sort {$b<=>$a} @aligned ) 
{
	$sum += $len ;
	print "$i $sum\n" ;
	++$i ;
}


