#!/bin/perl

# Compute the porprotion of covered area for a scaffold from the BAM depth file
# usage: a.pl xxx.fa yyy.depth

use strict ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;

my %length ;
my %covered ;
my %depth ;
my @cols ;
my $i ;
while ( <FP1> ) # Assume one-line fasta format
{
    chomp ;
    my $id = substr( $_, 1 ) ;
    my $seq = <FP1> ;
    chomp $seq ;
    @cols = split /N+/, $seq ;
    my $len = 0 ;
    for ( $i = 0 ; $i < scalar( @cols ) ; ++$i )
    {
            $len += length( $cols[$i] ) ;
    }
    # The length without gap
	$length{$id} = $len ;
 }
 
 while ( <FP2> )
 {
	chomp ; 
	@cols = split ;
	++$covered{ $cols[0] } ;
	$depth{ $cols[0] } += $cols[2] ;
 }
 
 foreach my $key (keys %length)
 {
	print "$key\t$length{$key}\t", $covered{$key} / $length{$key} * 100.0, "\t", $depth{$key} / $length{$key}, "\n" ;
 }
