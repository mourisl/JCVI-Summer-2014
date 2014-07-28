#!/bin/perl

# Compare the paired-end reads from apecca.bam and allpaths.bam
# usage: a.pl apecca.bam allpaths.bam

use strict ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;
my $line1 ;
my $prevLine1 = "" ;
my $line2 ;
my $prevLine2 = "" ;
my @reads1 ;
my @reads2 ;
my @cols1 ;
my @cols2 ;
my @cols ;
my $i ;

while ( <FP1> )
{
	# Collect a pair
	$line1 = $_ ;
	next if ( $line1 =~ /@/ ) ;
	if ( (split /\s/, $line1 )[0] eq ( split /\s/, $prevLine1 )[0] )
	{
		$prevLine1 = $line1 ;
		push @reads1, $line1 ;
	}
	else
	{
		while ( <FP2> )
		{
			$line2 = $_ ;
			next if ( $line2 =~ /@/ ) ;
			if ( (split/\s/, $line2)[0] eq ( split /\s/, $prevLine2 )[0] )
			{
				$prevLine2 = $line2 ;
				push @reads2, $line2 ;
			}
			else
			{
				last ;
			}
		}
		
		# Compare the reads
		my $tag1 = -1 ;
		my $tag2 = -1 ;
		
		for ( $i = 0 ; $i < scalar( @reads1 ) ; ++$i )
		{
		    if ( $reads1[$i] =~ /XP:Z:/ || $reads1[$i] =~ /\*/) 
			{
				$tag1 = -1 ;
				last ;
			}
			$tag1 = $i ;
		}
		for ( $i = 0 ; $i < scalar( @reads2 ) ; ++$i )
		{
			if ( $reads2[$i] =~ /XP:Z:/ || $reads2[$i] =~ /\*/ ) 
			{
				$tag2 = -1 ;
				last ;
			}
			$tag2 = $i ;
		}		
		
		if ( scalar( @reads1 ) != 2 || scalar( @reads2 ) != 2 )
		{
			$tag1 = -1 ;
			$tag2 = -1 ;
		}
		
		
		#print $reads1[$tag1], $reads2[$tag2] ;
		if ( $tag1 != -1 || $tag2 != -1 )
		{
			@cols1 = split /\s/, $reads1[$tag1] ;
			@cols2 = split /\s/, $reads2[$tag2] ;
			die "Unmatched read id\n" if ( $cols1[0] ne $cols2[0] ) ;
			my $result = 0 ;
			$result |= 1 if ( $cols1[6] eq "=" ) ;
			$result |= 2 if ( $cols2[6] eq "=" ) ; 
		
			#if ( $cols1[6] ne "*" && $cols1[6] ne "*" && 
			#	$result != 3 )
			#{
				print "# $result\n", $reads1[$tag1], $reads2[$tag2] ;
			#}
		}
		
		undef @reads1 ;
		push @reads1, $line1 ;
		undef @reads2 ;
		push @reads2, $line2 ;
		$prevLine2 = $line2 ;
	}
	$prevLine1 = $line1 ;
}
		while ( <FP2> )
		{
			$line2 = $_ ;
			next if ( $line2 =~ /@/ ) ;
			if ( (split/\s/, $line2)[0] eq ( split /\s/, $prevLine2 )[0] )
			{
				$prevLine2 = $line2 ;
				push @reads2, $line2 ;
			}
			else
			{
				last ;
			}
		}
	# Compare the reads
		my $tag1 = -1 ;
		my $tag2 = -1 ;
		
		for ( $i = 0 ; $i < scalar( @reads1 ) ; ++$i )
		{
		    if ( $reads1[$i] =~ /XP:Z:/ || $reads1[$i] =~ /\*/) 
			{
				$tag1 = -1 ;
				last ;
			}
			$tag1 = $i ;
		}
		for ( $i = 0 ; $i < scalar( @reads2 ) ; ++$i )
		{
			if ( $reads2[$i] =~ /XP:Z:/ || $reads2[$i] =~ /\*/ ) 
			{
				$tag2 = -1 ;
				last ;
			}
			$tag2 = $i ;
		}		
		
		if ( scalar( @reads1 ) != 2 || scalar( @reads2 ) != 2 )
		{
			$tag1 = -1 ;
			$tag2 = -1 ;
		}
		#print $reads1[$tag1], $reads2[$tag2] ;
		if ( $tag1 != -1 && $tag2 != -1 )
		{
			@cols1 = split /\s/, $reads1[$tag1] ;
			@cols2 = split /\s/, $reads2[$tag2] ;
			die "Unmatched read id\n" if ( $cols1[0] ne $cols2[0] ) ;
			# Core part:
			my $result = 0 ;
			$result |= 1 if ( $cols1[6] eq "=" ) ;
			$result |= 2 if ( $cols2[6] eq "=" ) ; 
		
			#if ( $cols1[6] ne "*" && $cols1[6] ne "*" && 
			#	$result != 3 )
			#{
				print "# $result\n", $reads1[$tag1], $reads2[$tag2] ;
			#}
		}