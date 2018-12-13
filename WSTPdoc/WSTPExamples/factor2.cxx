/* To run this program use the command-line below:
 *	Unix:           factor2 -linkname "math -wstp"
 *	Mac or Windows: factor2 -linkmode launch
 */

#include <stdio.h>
#include <stdlib.h>
#include "wstp.h"

static void init_and_openlink P((int argc, char* argv[]));
static void error P(( WSLINK lp));
static void read_and_print_expression P(( WSLINK lp));


WSENV ep = (WSENV)0;
WSLINK lp = (WSLINK)0;


int main( int argc, char* argv[])
{
	int n, pkt;

	init_and_openlink( argc, argv);

	printf( "Integer to factor: ");

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	scanf_s( "%d", &n);
#else
	scanf( "%d", &n);
#endif
	/* Send EvaluatePacket[ FactorInteger[n]]. */
	WSPutFunction( lp, "EvaluatePacket", 1L);
		WSPutFunction( lp, "FactorInteger", 1L);
			WSPutInteger( lp, n);
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


static void read_and_print_expression( WSLINK lp)
{
	const char *s;
	int n;
	int i, len;
	double r;
	static int indent;

	switch( WSGetNext( lp)) {
	case WSTKSYM:
		WSGetSymbol( lp, &s);
		printf( "%s ", s);
		WSReleaseSymbol( lp, s);
		break;
	case WSTKSTR:
		WSGetString( lp, &s);
		printf( "\"%s\" ", s);
		WSReleaseString( lp, s);
		break;
	case WSTKINT:
		WSGetInteger( lp, &n);
		printf( "%d ", n);
		break;
	case WSTKREAL:
		WSGetReal( lp, &r);
		printf( "%g ", r);
		break;
	case WSTKFUNC:
		indent += 3;
		printf( "\n %*.*s", indent, indent, "");
		if( WSGetArgCount( lp, &len) == 0){
			error( lp);
		}else{
			read_and_print_expression( lp);
			printf( "[");
			for( i = 1; i <= len; ++i){
				read_and_print_expression( lp);
				if( i != len) printf( ", ");
			}
			printf( "]");
		}
		indent -= 3;
		break;
	case WSTKERROR:
	default:
		error( lp);
	}
}


static void error( WSLINK lp)
{
	if( WSError( lp)) {
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

