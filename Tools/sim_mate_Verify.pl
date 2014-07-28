#!/bin/perl

# Verify the alignments of the sim_mate reads
# Usage: a.pl < in.sam

use strict ;
my @cols ;
my @lineCols ;
my $id ;
my $line ;
while (<>)
{
	$line = $_ ;
	@lineCols = split ;
	next if ( /XP:Z:/ ) ;
	
	if ( $lineCols[6] ne "=" && $lineCols[6] ne "*" )
	{
		print $_ ; 
	}
	
	next ;
	
	# Get the originated id and position
	$id = $lineCols[0] ;
	@cols = split /[#@]/, $id ;

	#print $cols[0], " ", $cols[1], "\n" ;
	my $scafId = "scaffold_".(split /_/, $cols[0])[2] ;
	my $fragStart = ( split /-/, $cols[1] )[0] ;
	my $fragEnd = ( split /-/, $cols[1] )[1] ; 
	
	# Get the aligned scafid and position
	my $alignId = $lineCols[2] ;
	my $start = $lineCols[3] ;
	my $cigar = $lineCols[5] ;
	my $softClip = 0 ;
	next if ( $cigar eq "*") ;
	my $AS = 2000 ;
	
	if ( $line =~ /AS:i:([0-9]+)?\s/ )
	{
		$AS = $1 ;
	}
	
	#next if ( $AS < 200 ) ;
	
	if ( $cigar =~ /^([0-9]+)S/ )
	{
		$softClip = $1 ;
		#print "$cigar $softClip\n" ;
	}
	$start -= $softClip ;
	if ( $alignId ne $scafId )
	{
		print $line ;
		next ;
	}
	
	if ( $lineCols[1] & 0x40 )
	{
		# Left mate
		if ( $start != $fragStart + 1 )
		{
			print $line ;
		}
	}
	else
	{	
		if ( $start + 1999 != $fragEnd )
		{
			print $line ;
		}
	}
}