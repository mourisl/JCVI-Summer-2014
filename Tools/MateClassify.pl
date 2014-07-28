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
		
		#print $reads1[$tag1], $reads2[$tag2] ;
		if ( defined $reads1[0]  && defined $reads2[0] )
		{
			my $result = 0 ;
			@cols1 = split /\s/, $reads1[0] ;
			@cols2 = split /\s/, $reads1[1] ;
			my $tmp = $cols1[0] ;
			#die "Unmatched read id\n" if ( $cols1[0] ne $cols2[0] ) ;
			$result |= 1 if ( $cols1[6] eq "=" && $cols2[6] eq "=" && 
								$cols1[5] ne "*" && $cols2[5] ne "*"
								&& scalar( @reads1 ) == 2 ) ;
			
			@cols1 = split /\s/, $reads2[0] ;
			@cols2 = split /\s/, $reads2[1] ;
			die "Unmatched read id\n" if ( $tmp ne $cols2[0] ) ;
		
			$result |= 2 if ( $cols1[6] eq "=" && $cols2[6] eq "=" && 
								$cols1[5] ne "*" && $cols2[5] ne "*"
								&& scalar( @reads2 ) == 2 ) ;
			
			@cols1 = split /\s/, $reads1[0] ;
			@cols2 = split /\s/, $reads2[0] ;
			print "# $result\n", 
					"$cols1[0]\t$cols1[2]\t$cols1[5]\t$cols1[6]\n", 
					"$cols2[0]\t$cols2[2]\t$cols2[5]\t$cols2[6]\n" ;
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
		my $result = 0 ;
		@cols1 = split /\s/, $reads1[0] ;
		@cols2 = split /\s/, $reads1[1] ;
		my $tmp = $cols1[0] ;
		#die "Unmatched read id\n" if ( $cols1[0] ne $cols2[0] ) ;
		$result |= 1 if ( $cols1[6] eq "=" && $cols2[6] eq "=" && 
								$cols1[5] ne "*" && $cols2[5] ne "*"
								&& scalar( @reads1 ) == 2 ) ;
		
		@cols1 = split /\s/, $reads2[0] ;
		@cols2 = split /\s/, $reads2[1] ;
		die "Unmatched read id\n" if ( $tmp ne $cols2[0] ) ;
	
		$result |= 2 if ( $cols1[6] eq "=" && $cols2[6] eq "=" && 
								$cols1[5] ne "*" && $cols2[5] ne "*"
								&& scalar( @reads2 ) == 2 ) ;
		
		@cols1 = split /\s/, $reads1[0] ;
		@cols2 = split /\s/, $reads2[0] ;
		print "# $result\n", 
				"$cols1[0]\t$cols1[2]\t$cols1[5]\t$cols1[6]\n", 
				"$cols2[0]\t$cols2[2]\t$cols2[5]\t$cols2[6]\n" ;