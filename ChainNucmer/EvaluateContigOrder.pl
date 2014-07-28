#!/bin/perl

# Choose the scaffold whose contig aligned well on both the reference and the scaffolds we want to test
my $help = "usage: a.pl xxx.contig.oneline.fa yyy.ref_xxx_contig.coords zzz.test_xxx_contig.coords 
						[-scafPortion sss | -verbose | -testGap test.oneline.fa | -refPortion rrr | -allowPartial]\n" ;

# output: the scaffold name and then the indices of contigs of interest in the order that the first number is smaller than the reverse order

use strict ;

use List::Util 'max' ;

die $help if ( scalar( @ARGV ) == 0 ) ;

open FP1, $ARGV[0] ;
open FP2, $ARGV[1] ;
open FP3, $ARGV[2] ;

my %contigLen ;
my %scafContigCnt ;
my %contigInRef ;
my %contigInTest ;
my %ignoreContig ;
my %scafLen ;
my %testScaf ;
my $scafPortion = 0 ;
my $verbose = 0 ;
my $testGap = 0 ;
my $refPortion = 0.75 ;
my $allowPartial = 0 ;
my $i ;
for ( $i = 3 ; $i < scalar( @ARGV ) ; ++$i )
{
	if ( $ARGV[$i] eq "-scafPortion" ) 
	{
		$scafPortion = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-verbose" )
	{
		$verbose = 1 ;
	}
	elsif ( $ARGV[$i] eq "-testGap" )
	{
		$testGap = 1 ;
		open FP4, $ARGV[$i + 1] ;
		while (<FP4>)
		{
			chomp ;
			my $id = (split /\s+/, substr( $_, 1 ))[0] ; 
			my $seq = <FP4> ;
			chomp $seq ;
			$testScaf{$id} = $seq ;
		}		
		close FP4 ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-refPortion" )
	{
		$refPortion = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-allowPartial" )
	{
		$allowPartial = 1 ;
	}
	else
	{
		die "Unknown argument $ARGV[$i]\n" ;
	}
}

sub Sign
{
	if ( $_[0] > 0 )
	{
		return "+" ;
	}
	elsif ( $_[0] < 0 )
	{
		return "-" ;
	}
	else
	{
		return "=" ;
	}
}

sub IsClose
{
	my $a = $_[0] ;
	my $b = $_[1] ;
	return 1 if ( $a > 1.2 * $b ) ;
	return 0 if ( $a >= 0.8 * $b && $a <= 1.2 * $b ) ;
	return -1 ;
}

sub ComputeLCS
{
	my $a = $_[0] ;
	my $b = $_[1] ;
	my $weight = $_[2] ;
	my @lcs0 ;
	my @lcs1 ;
	my $lena = scalar( @$a ) ;
	my $lenb = scalar( @$b ) ;
	my $i ;
	my $j ;
	
	$lcs0[0] = $lcs1[0] = 0 ;
	for ( $i = 0 ; $i < $lena ; ++$i )
	{
		for ( $j = 0 ; $j < $lenb ; ++$j )
		{
			if ( $a->[$i] eq $b->[$j] )
			{
				$lcs1[$j + 1] = $lcs0[$j] + $weight->{$a->[$i] } ;
			}
			else
			{
				$lcs1[$j + 1] = $lcs0[$j] ;
			}
			
			$lcs1[$j + 1] = $lcs0[$j + 1] if ( $lcs0[$j + 1] > $lcs1[$j + 1] ) ;
			$lcs1[$j + 1] = $lcs1[$j] if ( $lcs1[$j] > $lcs1[$j + 1] ) ;
		}
		@lcs0 = @lcs1 ;
	}
	return $lcs0[$lenb] ;
}

while (<FP1>)
{
	chomp ;
	my $id = (split /\s+/, substr( $_, 1 ))[0] ; 
	my $seq = <FP1> ;
	chomp $seq ;
	$contigLen{$id} = length( $seq ) ;
	
	++$scafContigCnt{ (split /\./, $id )[0] } ;
	$scafLen{ (split /\./, $id )[0] } += length( $seq ) ;
}

my $line ;
my @cols ;
#     3939     4038  |       25      124  |      100      100  |    99.00  | scaffold_202.56      scf620008933730
$line = <FP2> ;
$line = <FP2> ;
$line = <FP2> ;
$line = <FP2> ;
$line = <FP2> ;
while ( <FP2> )
{
	chomp ;
	$line = $_ ;
	@cols = split /\s+/, $line ;
	shift @cols if ( $cols[0] eq "" ) ;
	
	my $id = $cols[12] ;
	
	if ( defined $contigInRef{ $id } )
	{
		my $tmp = IsClose($cols[7], ${ $contigInRef{ $id } }[7] ) ;
		if (  $tmp == 1 )
		{
			@{ $contigInRef{ $id } } = @cols ;
			undef $ignoreContig{$id} ;
		}
		elsif ( $tmp == 0 )
		{
			@{ $contigInRef{ $id } } = @cols if ( $cols[7] > ${ $contigInRef{ $id } }[7] ) ;
			$ignoreContig{$id} = 1 ;
		}
	}
	else
	{
		@{ $contigInRef{ $id } } = @cols ;
	}
}

$line = <FP3> ;
$line = <FP3> ;
$line = <FP3> ;
$line = <FP3> ;
$line = <FP3> ;
while ( <FP3> )
{
	chomp ;
	$line = $_ ;
	@cols = split /\s+/, $line ;
	shift @cols if ( $cols[0] eq "" ) ;
	
	my $id = $cols[12] ;
	
	next if ( $ignoreContig{$id} == 1 ) ; # Masked by comparing with the reference.
	
	if ( defined $contigInTest{ $id } )
	{
		my $tmp = IsClose( $cols[7], ${ $contigInTest{ $id } }[7] ) ;
		if ( $tmp == 1 )
		{
			@{ $contigInTest{ $id } } = @cols ;
			undef $ignoreContig{$id} ;
		}
		elsif ( $tmp == 0 )
		{
			@{ $contigInTest{ $id } } = @cols if ( $cols[7] > ${ $contigInTest{ $id } }[7] ) ;
			$ignoreContig{$id} = 1 ;
		}
	}
	else
	{
		@{ $contigInTest{ $id } } = @cols ;
	}
}

# Evaluate each scaffolds.
my $sid ;
foreach $sid (keys %scafContigCnt)
{
	my $cnt = $scafContigCnt{ $sid } ;
	my $i ;
	my $id ;
	my $chosenContigLen ;
	my @choose ;
	my %startInRef ;
	my %startInTest ;
	my $refScaf = "" ;
	my $testScaf = "" ;
	#my $refDirection = 0 ;
	#my $testDirection = 0 ;
	my %refDirection ;
	my %testDirection ;
	
	# Decide use which scaffolds and chromosome
	# Based on where the longest scaffold aligned to
	my $max = 0 ;
	my $maxtag ;
	for ( $i = 0 ; $i < $cnt ; ++$i )
	{
		$id = $sid.".".$i ;
		my $len = $contigLen{$id} ;
		if ( $ignoreContig{$id} == 0 && $contigLen{ $id } > $max &&
			${ $contigInRef{$id} }[7] >= $refPortion * $len &&
			${ $contigInTest{ $id } }[7] >= 0.9 * $len )
		{
			$max = $len ;
			$maxtag = $id ;
		}
	}
	$refScaf = ${ $contigInRef{ $maxtag } }[11] ;
	$testScaf = ${ $contigInTest{ $maxtag } }[11] ;
	
	for ( $i = 0 ; $i < $cnt ; ++$i )
	{
		$id = $sid.".".$i ;
		next if ( $ignoreContig{ $id } == 1 ) ;
		my $len =  $contigLen{ $id } ;
		#next if ( $len < 1000 ) ;
		if ( ${ $contigInRef{$id} }[7] < $refPortion * $len ||
			${ $contigInTest{ $id } }[7] < 0.9 * $len )
		{
			if ( ${ $contigInRef{$id} }[7] >= 0.9 * $len 
				&& ${ $contigInTest{ $id } }[7] < 0.5 * $len 
				&& $verbose == 1 && $len >= 1000 ) 
			{
				print "WARNING: $id($len) misses in the other assembly.\n" ;
			}
			next ;
		}
		
		if ( $refScaf eq "" )
		{
			$refScaf = ${ $contigInRef{ $id } }[11] ;
		}
		elsif ( $refScaf ne ${ $contigInRef{ $id } }[11] )
		{
			#$i = -1 ;
			#last ;
			if ( $allowPartial == 1 )
			{
				next ;
			}
			else
			{
				$i = -1 ;
				last ;
			}
		}
		$refDirection{ $i } = Sign( ${ $contigInRef{$id} }[4] - ${ $contigInRef{$id} }[3] ) ;
		
		if ( $testScaf eq "" )
		{
			$testScaf = ${ $contigInTest{ $id } }[11] ;
		}
		elsif ( $testScaf ne ${ $contigInTest{ $id } }[11] )
		{
			#$i = -1 ;
			#last ;
			if ( $allowPartial == 1 )
			{
				next ;
			}
			else
			{
				$i = -1 ;
				last ;
			}
		}			
		$testDirection{$i} = Sign( ${ $contigInTest{$id} }[4] - ${ $contigInTest{$id} }[3] ) ;
		
		$startInRef{$i} = ${ $contigInRef{ $id} }[0] ;
		$startInTest{$i} = ${ $contigInTest{$id} }[0] ;
		
		$chosenContigLen += $len ;
		push @choose, $i ;
	}

	next if ( $i == -1 || scalar( @choose ) == 0 ) ;
	#print "WARNING: majority of the scaffold failed to align.\n" if ( $chosenContigLen < 0.1 * $scafLen{ $sid } ) ;
	next if ( $chosenContigLen < $scafPortion * $scafLen{ $sid } ) ;
	next if ( scalar( @choose ) < 2 ) ;
	print ( "\n$sid $refScaf $testScaf\n" ) ;
	my %contigIndLen ;
	for ( $i = 0 ; $i < scalar( @choose ) ; ++$i )
	{
		my $len = $contigLen{ $sid.".".($choose[$i]) } ;
		$contigIndLen{ $choose[$i] } = $len ;
		print $choose[$i],"(", $len, ") " ;
	}
	print "\n" ;
	
	my @refOrder ;
	my @refOrderPlus ;
	my @refOrderMinus ;
	sub byStartInRef
	{
		$startInRef{$a}<=>$startInRef{$b} ;
	}
	print " I: " ;
	foreach $i (sort {$startInRef{$a}<=>$startInRef{$b}} @choose) 
	{
		#print $startInRef{$i}, "\n" ;
		push @refOrder, $i ;
		if ( $refDirection{$i} eq "+" )
		{
			push @refOrderPlus, $i ;
		}
		else
		{
			push @refOrderMinus, $i ;
		}
		print "$i($refDirection{$i}) " ;
	}
	print "\n" ;
	
	my @testOrder ;
	my @testOrderPlus ;
	my @testOrderMinus ;
	sub byStartInTest
	{
		$startInTest{$a}<=>$startInTest{$b} ;
	}	
	print "II: " ;
	foreach $i (sort {$startInTest{$a}<=>$startInTest{$b}} @choose) 
	{
		push @testOrder, $i ;
		if ( $testDirection{$i} eq "+" )
		{
			push @testOrderPlus, $i ;
		}
		else
		{
			push @testOrderMinus, $i ;
		}		
		print "$i($testDirection{$i}) " ;
	}
	print "\n" ;
	
	my $lcsRef ;
	my $lcsTest ;
	my $lcsRefTest ;
	my $tmp ;
	my $refForward = 1 ;
	my $testForward = 1 ;
	my @tmpArray ;
	#next ;
	@refOrderMinus = reverse( @refOrderMinus ) ;
	@testOrderMinus = reverse( @testOrderMinus ) ;
	
	$lcsRef = ComputeLCS( \@choose, \@refOrderPlus, \%contigIndLen ) ;
	#@tmpArray = reverse(@refOrder ) ;
	$tmp = ComputeLCS( \@choose, \@refOrderMinus, \%contigIndLen ) ;
	if ( $tmp > $lcsRef )
	{
		$lcsRef = $tmp ;
		$refForward = 0 ;
	}

	$lcsTest = ComputeLCS( \@choose, \@testOrderPlus, \%contigIndLen ) ;
	#@tmpArray = reverse(@testOrder ) ;
	$tmp = ComputeLCS( \@choose, \@testOrderMinus, \%contigIndLen ) ;
	if ( $tmp > $lcsTest )
	{
		$lcsTest = $tmp  ;
		$testForward = 0 ;
	}
	
	@tmpArray = ( ComputeLCS( \@refOrderPlus, \@testOrderPlus, \%contigIndLen ),
					ComputeLCS( \@refOrderPlus, \@testOrderMinus, \%contigIndLen ),
					ComputeLCS( \@refOrderMinus, \@testOrderPlus, \%contigIndLen ), 
					ComputeLCS( \@refOrderMinus, \@testOrderMinus, \%contigIndLen ) 
				) ;
	$lcsRefTest = max( @tmpArray ) ;
	print "LCS: ", scalar(@choose)," ", $chosenContigLen, " $lcsRef $lcsTest $lcsRefTest\n" ;
	
	# Determine the events.
	if ( $verbose == 1 )
	{
		next if ( $lcsRef == $chosenContigLen && $lcsTest == $chosenContigLen ) ; # Everything agrees.
		if ( $lcsTest == $chosenContigLen )
		{
			my $len = 0 ;
			if ( $refForward == 1 && scalar( @refOrderMinus ) > 0 )
			{
				foreach $i ( @refOrderMinus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			elsif ( $refForward == 0 && scalar( @refOrderPlus ) > 0 )
			{
				foreach $i ( @refOrderPlus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			
			print "Event: Inversion in reference. Amount $len\n" if ( $len > 0 ) ;
		
			$len = 0 ;
			if ( $refForward == 1 )
			{
				foreach $i ( @refOrderPlus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			else
			{
				foreach $i ( @refOrderMinus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			#print "$len $lcsRef\n" ;
			print "Event: Translocation in reference. Amount ", $len - $lcsRef, "\n" if ( $len > $lcsRef ) ;
		}
		
		if ( $lcsRef == $chosenContigLen )
		{
			my $len = 0 ;
			if ( $testForward == 1 && scalar( @testOrderMinus ) > 0 )
			{
				foreach $i ( @testOrderMinus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			elsif ( $testForward == 0 && scalar( @testOrderPlus ) > 0 )
			{
				foreach $i ( @testOrderPlus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			
			print "Event: Inversion in reference. Amount $len\n" if ( $len > 0 ) ;
		
			$len = 0 ;
			if ( $testForward == 1 )
			{
				foreach $i ( @testOrderPlus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			else
			{
				foreach $i ( @testOrderMinus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			print "Event: Translocation in the other assembly. Amount ", $len - $lcsTest, "\n" if ( $len > $lcsTest ) ;
		}	
		
		if ( $lcsRefTest == $chosenContigLen )
		{
			my $len = 0 ;
			if ( $testForward == 1 && scalar( @testOrderMinus ) > 0 )
			{
				foreach $i ( @testOrderMinus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			elsif ( $testForward == 0 && scalar( @testOrderPlus ) > 0 )
			{
				foreach $i ( @testOrderPlus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			
			print "Event: Inversion in reference. Amount $len\n" if ( $len > 0 ) ;
		
			$len = 0 ;
			if ( $testForward == 1 )
			{
				foreach $i ( @testOrderPlus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}
			else
			{
				foreach $i ( @testOrderMinus )
				{
					$len += $contigIndLen{ $i } ;
				}
			}	
			if ( $len > $lcsTest ) 
			{
				print "Event: Translocation in the assembly. Amount ", $len - $lcsTest,  "\n" ;
			}
		}
		
		if ( $testGap == 1 )
		{
			for ( $i = 0 ; $i < scalar( @testOrder ) - 1 ; ++$i )
			{
				my $id = $sid.".".$testOrder[$i] ;
				my $id2 = $sid.".".$testOrder[$i + 1] ;
				my $from ;
				my $to ;
				my $containN = 0 ;
				$from = ${ $contigInTest{ $id } }[1] ;
				$to = ${ $contigInTest{ $id2 } }[0] ;
				if ( $from < $to )
				{
					#my $text =substr( $testScaf{ ${ $contigInTest{$id} }[11] }, $from - 1, $to - $from + 1 ) ;
					#print $text, "\n" ;
					$containN = 1 if ( substr( $testScaf{ ${ $contigInTest{$id} }[11] }, $from - 1, $to - $from + 1 ) =~ /N/ ) ;
				}
				print "Region: ${ $contigInTest{ $id } }[11] $from $to ", $to - $from + 1, " $containN\n" ;
			}
		}
	} # end if verbose
	print "\n" ;
}