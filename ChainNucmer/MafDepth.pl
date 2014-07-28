#!/bin/perl

my $usage = "a.pl xxx.oneline.fa yyy.maf" ;

die "$usage\n" if ( scalar( @ARGV ) == 0 ) ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my %scaf ;
my %origLength ;
my %coverage ;

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

while (<FP2>)
{
	chomp ;
	my $line = $_ ;
	my @cols ;
	my $i ;
	next if ( substr( $line, 0, 2 ) ne "a " ) ;
	$line = <FP2> ;
	chomp $line ;
	
	@cols = split /\s+/, $line ;
	
	my $array = \@{$coverage{ $cols[1] } };
	
	for ( $i = $cols[2] ; $i < $cols[2] + $cols[3] ; ++$i )
	{
		++( $array->[$i] ) ;
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

