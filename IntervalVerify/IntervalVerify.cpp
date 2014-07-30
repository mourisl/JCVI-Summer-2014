#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SCAFFOLD 100000

#define SAM_LINE_LEN 10000

char usage[] = "./a.out xxx.fa min_insert max_insert [-only_contig]< yyy.sam > output\n"
			   "The input yyy.sam should be sorted by read id" ;
char buffer[5000000] ;

bool onlyContig = false ;

struct _tree // segment tree
{
	unsigned int cov ;
	int start, end ;
	struct _tree *left, *right ;
} ;

struct _scaffold 
{
	char name[100] ;
	char *seq ;
	int len ;
	int *coverage ;
	struct _tree *tree ;
} ;

struct _scaffold scaffolds[MAX_SCAFFOLD] ; 


struct _mate
{
	char name[SAM_LINE_LEN] ;
	int readCnt ;
	char reads[20][SAM_LINE_LEN] ;
} ;

struct _tree *BuildTree( int s, int e )
{
	struct _tree *ret = ( struct _tree * )malloc( sizeof( *ret ) ) ;
	ret->start = s ;
	ret->end = e ;
	ret->cov = 0 ;
	
	if ( e == s )
	{
		ret->left = NULL ;
		ret->right = NULL ;
		return ret ;
	}
	
	ret->left = BuildTree( s, ( e + s ) / 2 ) ;
	ret->right = BuildTree( ( e + s ) / 2 + 1, e ) ;
	
	return ret ;
}

void InsertInterval( struct _tree *node, int s, int e ) 
{
	if ( s <= node->start && e >= node->end )
	{
		++node->cov ;
		return ;
	}
	if ( s <= ( node->start + node->end ) / 2 )
		InsertInterval( node->left, s, e ) ;
	if ( e > ( node->start + node->end ) / 2 )
		InsertInterval( node->right, s, e ) ;
}

int PositionCoverage( struct _tree *node, int pos )
{
	if ( node == NULL )
		return 0 ;
	int ret = node->cov ;
	if ( pos < ( node->start + node->end ) / 2 )
		ret += PositionCoverage( node->left, pos ) ;
	else
		ret += PositionCoverage( node->right, pos ) ;
	return ret ;
}

int CompScafId( const void *p1, const void *p2 )
{
	return strcmp( (( struct _scaffold *)p1)->name, ( (struct _scaffold *)p2)->name ) ;
}

void Swap( int &a, int &b )
{
	int tmp = b ;
	b = a ;
	a = tmp ;
}

void ProcessMate( const struct _mate &mate, int minInsert, int maxInsert, struct _scaffold scaffolds[], int scafCnt ) 
{
	if ( mate.readCnt > 2 ) // Ignore chimeric alignments
		return ;
	int tag1, tag2, pos1, pos2 ;
	char scafId1[SAM_LINE_LEN], scafId2[SAM_LINE_LEN] ;
	int l, r, m ;
	
	sscanf( mate.reads[0], "%s %d %s %d", buffer, &tag1, scafId1, &pos1 ) ;
	sscanf( mate.reads[1], "%s %d %s %d", buffer, &tag2, scafId2, &pos2 ) ;
	
	if ( strcmp( scafId1, scafId2 ) )
		return ;
	
	l = 0 ;
	r = scafCnt - 1 ;
	
	while ( l <= r )
	{
		m = ( l + r ) / 2 ;
		int tmp = strcmp( scafId1, scaffolds[m].name ) ; 
		if ( tmp == 0 )
			break ;
		else if ( tmp < 0 )
			r = m - 1 ;
		else
			l = m + 1 ;
	}
	
	if ( l > r )
	{
		//printf( "Unknown scaffolid name %s found in SAM\n", scafId1 ) ;
		//exit( 1 ) ;
		return ;
	}

	// Test whether it is Innie:
	if ( ( tag1 & 0x10 ) == ( tag2 & 0x10 ) )
		return ;
	
	if ( ( ( tag1 & 0x10 ) != 0 ) == ( pos1 < pos2 ))
		return ;
	//if ( ( ( tag1 & 0x10 ) != 0 ) == ( pos1 > pos2 ) )
	//	return ;
	
	// Test the insert size
	int insert, start, end ;	
	sscanf( mate.reads[0], "%s%s%s%s%s%s%s%s%d", buffer, buffer, buffer, buffer, buffer, 
								  buffer, buffer, buffer, &insert ) ;
	if ( insert < 0 )
		insert = -insert ;
	start = ( pos1 < pos2 ? pos1 : pos2 ) ;
	end = start + insert - 1 ;
	
	// Hard-coded clip
	start += 100 ;
	end -= 100 ;
	if ( end < start )
		return ;

	if ( insert >= minInsert && insert <= maxInsert )
	{
		//if ( !strcmp( mate.name, "HISEQ2:1018:C3YMPACXX:8:1101:8327:4993" ) )
		//	printf( "%s\n%s\n%d %d\n", mate.reads[0], mate.reads[1], start, end ) ;	
		InsertInterval( scaffolds[m].tree, start, end ) ;
	}
}

void VerifyScaffold( struct _scaffold &scaf )
{
	char *v = ( char * )malloc( sizeof( char ) * ( scaf.len + 10 ) ) ;
	int i ;
	for ( i = 0 ; i < scaf.len ; ++i )
	{
		int tmp = PositionCoverage( scaf.tree, i + 1 ) ;
		//if ( i == 290000 )
		//	printf( "%d\n", tmp ) ;
		if ( tmp > 0 )
			v[i] = '1' ;
		else
			v[i] = '0' ;
	}
	
	if ( onlyContig )
	{
		char pattern[21] ;
		char *p ;
		for ( i = 0 ; i < 20 ; ++i )
			pattern[i] = 'N' ;
		pattern[i] = '\0' ;
		p = scaf.seq ;
		while ( 1 )  
		{
			p = strstr( p, pattern ) ;
			if ( p == NULL )
				break ;
			while ( *p == 'N' )
			{
				v[p - scaf.seq] = '0' ;
				++p ;
			}
		}		
	}
	
	v[i] = '\0' ;
	printf( ">%s\n%s\n", scaf.name, v ) ;
	free( v ) ;
}

int CompareId( char *a, char *b )
{
	int la = strlen( a ) ;
	int lb = strlen( b ) ;
	if ( la <= 2 || lb <=2 || la != lb )
		return strcmp( a, b ) ;
	if ( a[la - 2] == '#' && b[lb - 2] == '#' )
	{
		a[la - 1] = '#' ;
		b[lb - 1] = '#' ;
	}
	return strcmp( a, b ) ;
}

int main( int argc, char *argv[] ) 
{
	int i, j, k ;
	int minInsert, maxInsert ;
	FILE *fpFa ;
	int scafCnt = 0 ;
	int scafLen = 0 ;
	char line[SAM_LINE_LEN] = "" ;
	
	if ( argc < 4 )
	{
		printf( "%s\n", usage ) ;
		return 1 ;
	}
	
	onlyContig = false ;
	for ( i = 0 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-only_contig" ) ) 
			onlyContig = true ;
	}
	
	// Get the size of scaffolds
	fpFa = fopen( argv[1], "r" ) ;
	while ( fgets( buffer, sizeof( buffer ), fpFa ) != NULL )
	{
		if ( buffer[0] == '>' )
		{
			if ( scafCnt > 0 )
				scaffolds[ scafCnt - 1 ].len = scafLen ;		
			
			sscanf( buffer + 1, "%s", scaffolds[ scafCnt ].name ) ;
			++scafCnt ;
			scafLen = 0 ;
		}
		else
		{
			int len = strlen( buffer ) ;
			if ( buffer[len - 1] == '\n' )
				--len ;
			scafLen += len ;
		}
	}
	if ( scafCnt > 0 )
		scaffolds[ scafCnt - 1 ].len = scafLen ;
		
	fclose( fpFa ) ;
	
	// Read in the sequence of scaffolds;
	fpFa = fopen( argv[1], "r" ) ;
	i = 0 ;
	while ( fgets( buffer, sizeof( buffer ), fpFa ) != NULL )
	{
		if ( buffer[0] == '>' )
		{
			scaffolds[i].seq = ( char * )malloc( sizeof( char ) * ( scaffolds[i].len + 10 ) ) ;
			scaffolds[i].seq[0] = '\0' ;
			++i ;
		}
		else
		{
			int len = strlen( buffer ) ;
			if ( buffer[len - 1] == '\n' )
				buffer[len - 1] = '\0' ;
			strcat( scaffolds[i - 1].seq, buffer ) ;
		}
	}
	
	// Sort the scaffolds by their names.
	qsort( scaffolds, scafCnt, sizeof( scaffolds[0] ), CompScafId ) ;
	
	// Build the segment tree for each scaffold
	for ( i = 0 ; i < scafCnt ; ++i )
	{
		scaffolds[i].tree = BuildTree( 1, scaffolds[i].len ) ;
	}
	
	minInsert = atoi( argv[2] ) ;
	maxInsert = atoi( argv[3] ) ;

	struct _mate mate ;
	mate.readCnt = 0 ;

	while ( fgets( line, sizeof( line ), stdin ) != NULL )
	{
		if ( line[0] == '@' )
			continue ;
			
		char id[SAM_LINE_LEN] ;

		sscanf( line, "%s", id ) ;
	
		if ( mate.readCnt > 0 && CompareId( id, mate.name ) )
		{
			ProcessMate( mate, minInsert, maxInsert, scaffolds, scafCnt ) ;
			mate.readCnt = 1 ;
			strcpy( mate.name, id ) ;
			strcpy( mate.reads[0], line ) ;
		}
		else
		{
			if ( mate.readCnt == 0 )
				strcpy( mate.name, id ) ;
			if ( mate.readCnt < 20 )
			{
				strcpy( mate.reads[ mate.readCnt ], line ) ;
				++mate.readCnt ;
			}
		}
	}
	ProcessMate( mate, minInsert, maxInsert, scaffolds, scafCnt ) ;
	
	for ( i = 0 ; i < scafCnt ; ++i )
		VerifyScaffold( scaffolds[i] ) ;
	return 0 ;
}
