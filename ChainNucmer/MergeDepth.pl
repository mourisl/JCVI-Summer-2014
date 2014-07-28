#!/bin/perl
# Merge a bunch of depths file. Also add the N positions

# usage: a.pl xxx.fa depth_file_list

use strict ;

my %origLength ;
my %scaf ;
my %coverage ;
my $i ;
open FP1, (shift @ARGV ) ;

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
close FP1 ;

while (<>)
{
	my @cols ;
	chomp ;
	@cols = split ;
	if ( $cols[2] =~ /[0-9]/ )
	{
		${$coverage{ $cols[0] }}[$cols[1]] += $cols[2] ;
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