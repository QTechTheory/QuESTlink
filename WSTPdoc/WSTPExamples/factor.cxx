/* To run this program use the command-line below:
 *	Unix:           factor -linkname "math -wstp"
 *	Mac or Windows: factor -linkmode launch
 */


#include <stdio.h>
#include <stdlib.h>

#include "wstp.h"

static void init_and_openlink( int argc, char* argv[]);
static void error( WSLINK lp);


WSENV ep = (WSENV)0;
WSLINK lp = (WSLINK)0;


int main(int argc, char* argv[])
{
	int pkt, n, prime, expt;
	int len, lenp, k;

	init_and_openlink( argc, argv);

	printf( "Integer to factor: ");

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	scanf_s( "%d", &n);
#else
	scanf( "%d", &n);
#endif

	WSPutFunction( lp, "EvaluatePacket", 1L);
		WSPutFunction( lp, "FactorInteger", 1L);
			WSPutInteger( lp, n);
	WSEndPacket( lp);

	while( (pkt = WSNextPacket( lp), pkt) && pkt != RETURNPKT) {
		WSNewPacket( lp);
		if (WSError( lp)) error( lp);
	}

	if ( ! WSTestHead( lp, "List", &len)) error(lp);
	for (k = 1; k <= len; k++) {
		if (WSTestHead( lp, "List", &lenp)
		&&  lenp == 2
		&&  WSGetInteger( lp, &prime)
		&&  WSGetInteger( lp, &expt)
		){
			printf( "%d ^ %d\n", prime, expt);
		}else{
			error( lp);
		}
	}

	WSPutFunction( lp, "Exit", 0);

	return 0;
}


static void error( WSLINK lp)
{
	if( WSError( lp)){
		fprintf( stderr, "Error detected by WSTP: %s.\n",
			WSErrorMessage(lp));
	}else{
		fprintf( stderr, "Error detected by this program.\n");
	}
	exit(3);
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
