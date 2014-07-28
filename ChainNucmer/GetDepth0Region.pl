#!/bin/perl

# Output the region with depth 0
# usage: a.pl xxx.oneline.fa yyy.depth

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my %origLength ;
my %scaf ;
while ( <FP1> ) # Assume one-line fasta format
{
	chomp ;
	my $id = (split /\s+/, substr( $_, 1 ))[0] ; 
	my $seq = <FP1> ;
	chomp $seq ;

	#$origLength{$id} = length( $seq ) ;
	$origLength{$id} = length( $seq ) ;
}
#print "phase 1 finished\n" ;
my $line = "" ;
my $prevLine = "" ;
my @cols1 ;
my @cols2 ;
while ( <FP2> )
{
	chomp ;
	my $len ;
	$line = $_ ;
	@cols1 = split /\s+/, $line ;
	next if ( $cols1[2] eq "0" ) ;
	#@cols2 = split /\s+/, $prevLine ;
	if ( $cols1[0] ne $cols2[0] )
	{
		if ( defined $cols2[0] && $cols2[1] < $origLength{ $cols2[0] } )
		{
			$len = $origLength{ $cols2[0] } - $cols2[1] ;
			my $type = "ES" ;
			$type = "WS" if ( $cols2[1] == 0 ) ;
			print "$cols2[0] ".($cols2[1]+1)." $origLength{ $cols2[0] } $len $type\n" ;
		}
		if ( defined $cols1[0] && $cols1[1] > 1 )
		{
			print "$cols1[0] 1 ".($cols1[1] - 1)." ".( $cols1[1] - 1)." SS\n" ;
		}
	}
	else
	{
		if ( $cols1[1] > $cols2[1] + 1 )
		{
			$len = $cols1[1] - $cols2[1] - 1 ;
			my $type = "MC" ;
			if ( $cols1[2] eq "N" )
			{
				$type = "EC" ;
			}
			if ( $cols2[2] eq "N" )
			{
				if ( $type eq "EC" )
				{
					$type = "WC" ;
				}
				else
				{
					$type = "SC" ;
				}
			}
			print "$cols1[0] ".($cols2[1] + 1)." ".($cols1[1] - 1)." $len $type\n" ;
		}
	}
	@cols2 = @cols1 ;
}

#@cols2 = split /\s+/, $prevLine ;
if ( defined $cols2[0] && $cols2[1] < $origLength{ $cols2[0] } )
{
	$len = $origLength{ $cols2[0] } - $cols2[1] ;
	my $type = "ES" ;
	$type = "WS" if ( $cols2[1] == 0 ) ;
	print "$cols2[0] ".($cols2[1]+1)." $origLength{ $cols2[0] } $len $type\n" ;
}