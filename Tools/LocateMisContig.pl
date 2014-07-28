#!/bin/perl

# Use the chimerical alignments to locate possible misassembled contigs
# usage: perl a.pl < xxx.sam

#HISEQ2:1018:C3YMPACXX:8:1101:1807:1954  161   scf620008936568   374833  60  56S44M     =   374855  122  
#AGATATGTTCAAATAAGGCCAAACATGCCTTAAAGCTAGATGTGTATAAGAGACAGGTATTAGAATGACCCTTGTCGTGAAAATTAGATTTTAGCCTATA  
#CCCFFFFFHHHGHJJJJJJJJIIJJJJJJJJJJJJJJJJJJJDHIIIIJJJJJJJJJFHIIJIJIIJJJJJJJJIJJHH;?CE7?>;@ADCE;35@C:A;    
#NM:i:0MD:Z:44 AS:i:44 XS:i:34 SA:Z:scf620008936568,383905,+,37M63S,60,0;


use strict ;
my %chimericInfo ;
my $line ;
my @cols ;
my @otherAlignment ;
my %align ;
my %alignEnd ;
my $i ;
my $text ;
my $SA ;
my $key ;
my @tmpCols ;
my $span ;
my $strand ;

# Return the start offset w.r.t the read on the reference direction 
# and the spanning of this alignment on the reference
sub ParseCigar
{
	my $cigar = $_[0] ;
	my $strand = $_[1] ;
	my @op ;
	my @num ;
	my $start ;
	my $span ;
	my $i ;
	@op = split /[0-9]+/, $cigar ;
	shift @op ;
	@num = split /[A-Za-z]/, $cigar ;
	#print @op, "\n" ;
	#print @num, "\n" ;
	
	if ( $strand eq "-" )
	{
		@op = reverse( @op ) ;
		@num = reverse( @num ) ;
	}
	
	$start = 0 ;
	$start = $num[0] if ( $op[1] eq "S" ) ;
	$i = 0 ;
	$i = 1 if ( $op[0] eq "S" ) ;
	$span = 0 ;
	for ( ; $i < scalar( @op ) ; ++$i )
	{
		last if ( $op[$i] eq "S" ) ;
		$span += $num[$i] ;
		$span -= $num[$i] if ( $op[$i] eq "I" ) ;
	}
	#print $num[0], " ", , "\n" ;
	return ($start, $span) ;
}

while ( <STDIN> )
{
	next if ( /^@/ || !( /SA:Z:/ ) ) ;	
	chomp ;
	$line = $_ ;
	@cols = split /\s+/, $line ;
	next if ( ( $cols[1] & 0x800 ) != 0 ) ;
	if ( ( $cols[1] & 0x10 ) != 0 )
	{
		$strand = '-' ;
	}
	else
	{
		$strand = '+' ;
	}
	
	# For a read, $key is the start position of alignment in the read
	($key, $span) = ParseCigar( $cols[5], $strand ) ;
	$text = "$cols[2]\t$cols[3]" ;
	undef %align ;
	undef %alignEnd ;
	#print "0: $key\n" ;
	$align{ $key } = $text ;
	$text = "$cols[2]\t".($cols[3] + $span - 1) ;
	$alignEnd{$key} = $text ;
	
	if ( $line =~ /SA:Z:(.+)?/ )
	{
		$SA = ( split/\s/, $1 )[0];
		@otherAlignment = split /;/, $SA ;
		for ( $i = 0 ; $i < scalar( @otherAlignment ) ; ++$i )
		{
			@cols = split /,/, $otherAlignment[$i] ;
			($key, $span) = ParseCigar( $cols[3], $cols[2] ) ;
			$text = "$cols[0]\t$cols[1]" ;
			#print "1: $key\n" ;
			#print $cols[0], "\t", $text, "\n" ;
			$align{ $key } = $text ;
			$text = "$cols[0]\t".($cols[1] + $span - 1) ;	
			$alignEnd{ $key } = $text ;
		}
	}
	else
	{
		next ;
	}
	
	# Update the chimeric part
	my $prevKey = -1 ;
	foreach $key (sort {$a<=>$b} (keys %align) )
	{
		if ( $prevKey == -1 ) 
		{
			$prevKey = $key ;
			next ;
		}
		if ( $alignEnd{ $prevKey } < $align{$key} )
		{
			$text = $alignEnd{$prevKey}."\t".$align{$key} ;
		}
		else
		{
			$text = $alignEnd{$prevKey}."\t".$align{$key} ;
		}
		#print $text, "\n" ;
		if ( defined $chimericInfo{$text} )
		{
			++$chimericInfo{$text} ;
		}
		else
		{
			$chimericInfo{$text} = 1 ;
		}
		$prevKey = $key ;
	}
}
# Output the result
foreach $key (keys %chimericInfo)
{
	print "$key\t$chimericInfo{$key}\n" ;
}