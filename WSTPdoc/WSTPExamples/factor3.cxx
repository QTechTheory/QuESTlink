/* To run this program use the command-line below:
 *	Unix:           factor3 -linkname "math -wstp"
 *	Mac or Windows: factor3 -linkmode launch
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wstp.h"

static void init_and_openlink P(( int argc, char* argv[]));
static void error P(( WSLINK lp));
static int   read_and_print_expression P(( WSLINK lp));
static int   read_and_print_atom P(( WSLINK lp, int tag));
static int   read_and_print_function P(( WSLINK lp));


WSENV ep = (WSENV)0;
WSLINK lp = (WSLINK)0;


int main( int argc, char* argv[])
{
	char buf[BUFSIZ];
	int pkt;

	init_and_openlink( argc, argv);

	printf( "Integer to factor: ");

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	scanf_s( "%s", buf, BUFSIZ);
#else
	scanf( "%s", buf);
#endif

	/* Send EvaluatePacket[ FactorInteger[n]]. */
	WSPutFunction( lp, "EvaluatePacket", 1L);
		WSPutFunction( lp, "FactorInteger", 1L);
			WSPutNext( lp, WSTKINT);
			WSPutByteString( lp, (unsigned char *)buf, (long)strlen(buf));
	WSEndPacket( lp);


	/* skip any packets before the first ReturnPacket */
	while( (pkt = WSNextPacket( lp), pkt) && pkt != RETURNPKT) {
		WSNewPacket( lp);
		if( WSError( lp)) error( lp);
	}

	read_and_print_expression( lp);
	printf( "\n");

	WSPutFunction( lp, "Exit", 0L);

	return 0;
}


static int read_and_print_expression( WSLINK lp)
{
	int tag;

	switch (tag = WSGetNext( lp)) {
	case WSTKSYM:
	case WSTKSTR:
	case WSTKINT:
	case WSTKREAL:
		return read_and_print_atom( lp, tag);
	case WSTKFUNC:
		return (read_and_print_function( lp));
	case WSTKERROR:
	default:
		return 0;
	}
}


static int   read_and_print_atom( WSLINK lp, int tag)
{
	const char *s;
	if( tag == WSTKSTR) putchar( '"');
	if( WSGetString( lp, &s)){
		printf( "%s", s);
		WSReleaseString( lp, s);
	}
	if( tag == WSTKSTR) putchar( '"');
	putchar( ' ');
	return WSError( lp) == WSEOK;
}


static int read_and_print_function( WSLINK lp)
{
	int  len, i;
	static int indent;

	if( ! WSGetArgCount( lp, &len)) return 0;

	indent += 3;
	printf( "\n%*.*s", indent, indent, "");

	if( read_and_print_expression( lp) == 0) return 0;
	printf( "[");

	for( i = 1; i <= len; ++i) {
		if( read_and_print_expression( lp) == 0) return 0;
		if( i < len) printf( ", ");
	}
	printf( "]");
	indent -= 3;

	return 1;
}


static void error( WSLINK lp)
{
	if (WSError( lp)) {
		fprintf( stderr, "Error detected by WSTP: %s.\n",
		WSErrorMessage( lp));
	}else{
		fprintf( stderr, "Error detected by this program.\n");
	}
	exit( 1);
}


static void deinit( void)
{
	if( ep) WSDeinitialize( ep);
}


static void closelink( void)
{
	if( lp) WSClose( lp);
}


static void init_and_openlink( int argc, char* argv[])
{
	int err;

	ep =  WSInitialize( (WSParametersPointer)0);
	if( ep == (WSENV)0) exit(1);
	atexit( deinit);

	lp = WSOpenArgv( ep, argv, argv + argc, &err);
	if(lp == (WSLINK)0) exit(2);
	atexit( closelink);	
}
