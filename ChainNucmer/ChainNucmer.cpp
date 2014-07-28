// Chain the local alignment of Nucmer 
/**
 The format of the input file should be:
 [S: name of the contig or scaffold] [Length of S] [n: # of local aligments for a contig or scaffold we want to chain]
  n lines: each line start and end  position of S, name of scaffold that this part of S aligned to, start and end of the destination
 **/
 // Usage: ./a.out [-noself | -crossScaffold | -onlyself |-showChain | -overlap xxx | distance xxx]< in

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FRAGMENT 65536
#define MAX_ID_LEN 64

int dp[MAX_FRAGMENT] ;
int dpChoose[MAX_FRAGMENT] ;
int tmpChoose[MAX_FRAGMENT], choose[MAX_FRAGMENT] ;
bool noSelfAlignment = false ;
bool onlySelfAlignment = false ;
bool showChain = false ;
int overlap = 0 ; // Allow overlap of two fragments
int distance = -1 ; // The difference of adjacent in two genomes should not be larger than this.

struct _fragment
{
	int start, end ;
	char destId[MAX_ID_LEN] ;
	int dstart, dend ;
} ;

int Comp( const void *p1, const void *p2 )
{
	struct _fragment *f1 = (struct _fragment *)p1 ;
	struct _fragment *f2 = (struct _fragment *)p2 ;
	int tmp = strcmp( f1->destId, f2->destId ) ;
	if ( tmp != 0 )
		return tmp ;
	else 
		return f1->start - f2->start ;
}

bool SameSign( int a, int b )
{
	// Avoid overflow
	if ( a > 0 && b > 0 )
		return true ;
	if ( a < 0 && b < 0 )
		return true ;
	return false ;
}

int Dp( int tag, struct _fragment *fragments, int n, char *id )
{
	if ( dp[tag] != 0 )
		return dp[tag] ;
	int i ;
	
	if ( ( noSelfAlignment || onlySelfAlignment ) && !strcmp( id, fragments[tag].destId ) 
		&& ( fragments[tag].start ==  fragments[tag].dstart || fragments[tag].end == fragments[tag].dend ) )
	{
		dpChoose[tag] = -1 ;
		return dp[tag] = -1 ;
	}
	
	int ret = fragments[tag].end - fragments[tag].start + 1 ;
	int max = 0 ;
	int maxtag = -1 ;
	for ( i = tag + 1 ; i < n ; ++i )
	{
		if ( fragments[i].start > fragments[tag].end - overlap && fragments[i].start > fragments[tag].start
			&& SameSign( fragments[i].dend - fragments[i].dstart, fragments[tag].dend - fragments[tag].dstart ) )
		{
			int d1 = fragments[i].start - fragments[tag].end ;
			int d2 ;
			// The coordinate on the destination should also be monotonic
			if ( fragments[tag].dend - fragments[tag].dstart > 0 )
			{
				if ( fragments[i].dstart <= fragments[tag].dend - overlap )
					continue ;
				d2 = fragments[i].dstart - fragments[tag].dend ;
			}
			else
			{
				if ( fragments[i].dstart >= fragments[tag].dend + overlap )
					continue ;
				d2 = fragments[tag].dend - fragments[i].dstart ;
			}
			
			if ( distance != -1 && ( d1 - d2 > distance || d2 - d1 > distance ) ) 
				continue ;
			
			int tmp = Dp( i, fragments, n, id ) ;
			if ( tmp > max )
			{
				max = tmp ;
				maxtag = i ;
			}
		}
	}
	
	dpChoose[tag] = maxtag ;
	return dp[tag] = ret + max ;
}

// Return the total length of the chaining
int Chain( struct _fragment *fragments, int n, char *id, int *choose )
{
	//printf( "%d\n", n ) ;
	if ( onlySelfAlignment && strcmp( id, fragments[0].destId ) )
		return -1 ;
	memset( dp, 0, sizeof( dp ) ) ;
	memset( dpChoose, -1, sizeof( dpChoose ) ) ;
	int max = -1 ;
	int maxtag = -1 ;
	int i ;
	for ( i = 0 ; i < n ; ++i )
		Dp( i, fragments, n, id ) ;
	for ( i = 0 ; i < n ; ++i )
		if ( dp[i] > max )
		{
			max = dp[i] ;
			maxtag = i ;
		}
		
	int p = maxtag ;
	choose[0] = p ;
	for ( i = 1 ; p != -1 ; ++i )
	{
		p = dpChoose[p] ;
		choose[i] = p ;
	}
	return max ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	char id[MAX_ID_LEN], destId[MAX_ID_LEN] ;
	int lenS, n ;
	struct _fragment *fragments ;
	int tmp, max, maxtag ;

	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-noself" ) )
			noSelfAlignment = true ;
		else if ( !strcmp( argv[i], "-onlyself" ) )
			onlySelfAlignment = true ;
		else if ( !strcmp( argv[i], "-showChain" ) )
			showChain = true ;
		else if ( !strcmp( argv[i], "-overlap" ) )
		{
			overlap = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( argv[i], "-distance" ) )
		{
			distance = atoi( argv[i + 1] ) ;
			++i ;
		}
		else
		{
			printf( "Unknown parameter %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}
	
	while ( scanf( "%s %d %d", id, &lenS, &n ) != EOF )
	{
		//printf( "++++ %d\n", n ) ;
		fragments = (struct _fragment *)malloc( sizeof( *fragments ) * n ) ;
		for ( i = 0 ; i < n ; ++i )
		{
			struct _fragment &f = fragments[i] ;
			scanf( "%d %d %s %d %d", &f.start, &f.end, f.destId, &f.dstart, &f.dend ) ;
		}
		qsort( fragments, n, sizeof( *fragments ), Comp ) ;
		
		int tag = 0 ;
		max = 0 ;
		maxtag = -1 ;
		for ( i = 1 ; i < n ; ++i )
		{
			if ( strcmp( fragments[i].destId, fragments[tag].destId ) )
			{
				/*if ( noSelfAlignment && !strcmp( fragments[tag].destId, id ) )
				{
					tag = i ;
					continue ;
				}*/
				tmp = Chain( fragments + tag, i - tag, id, tmpChoose ) ;
				if ( tmp > max )
				{
					max = tmp ;
					maxtag = tag ;
					
					int c ;
					for ( c = 0 ; tmpChoose[c] != -1 ; ++c )
						choose[c] = tmpChoose[c] ;
					choose[c] = -1 ;
				}
				tag = i ;
			}
		}
		//if ( !( noSelfAlignment && !strcmp( fragments[tag].destId, id ) ) )
		//{
			tmp = Chain( fragments + tag, n - tag, id, tmpChoose ) ;
			if ( tmp > max )
			{
				max = tmp ;
				maxtag = tag ;
				
				int c ;
				for ( c = 0 ; tmpChoose[c] != -1 ; ++c )
					choose[c] = tmpChoose[c] ;
				choose[c] = -1 ;				
			}
		//}
		if ( maxtag == -1 )
			strcpy( destId, "-") ;
		else
			strcpy( destId, fragments[maxtag].destId ) ;
			
		printf( "%s %d %d %lf %s", id, max, lenS, (double)max / lenS * 100.0, destId ) ;
		
		if ( showChain )
		{
			for ( i = 0 ; choose[i] != -1 && maxtag != -1 ; ++i )
				;
			printf( " %d", i ) ;
			for ( i = 0 ; choose[i] != -1 && maxtag != -1 ; ++i )
				printf( " [(%d %d) (%d %d)]", fragments[ choose[i] + maxtag ].start, fragments[ choose[i] + maxtag ].end,
					fragments[ choose[i] + maxtag ].dstart, fragments[ choose[i] + maxtag ].dend ) ;		
			printf( "\n" ) ;
		}
		else
			printf( "\n" ) ;
		free( fragments ) ;
	}
	return 0 ;
}