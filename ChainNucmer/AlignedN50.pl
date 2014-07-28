#!/bin/perl

my @len ;
my @cols ;
my $sum = 0 ;
while (<stdin>)
{
	chomp ;
	@cols = split ;
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
		push @len, $end - $start + 1 ;
		$sum += ( $end - $start + 1) ;
	}
}

my $psum = 0 ;
my $i = 0 ;
$sum = 107992474 ;
print "N stats with respect to $sum\n" ;
foreach my $l (sort {$b <=> $a} @len )
{
	$psum += $l ;
	while( $psum >= $sum * $i / 10 )
	{
		print "N".($i * 10)." $l\n" ;
		++$i ;
	}
}