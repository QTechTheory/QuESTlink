/* To run this program use the command-line below:
 *	Unix:           quotient -linkname "math -wstp"
 *	Mac or Windows: quotient -linkmode launch
 */

#include <stdio.h>
#include <stdlib.h>
#include "wstp.h"

static void init_and_openlink( int argc, char* argv[]);



WSENV ep = (WSENV)0;
WSLINK lp = (WSLINK)0;


int main( int argc, char* argv[])
{
	int pkt, m, n, q;

	init_and_openlink( argc, argv);

	printf( "Two integers m/n: ");

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	if( scanf_s( "%d/%d", &m, &n) != 2 && scanf_s( "%d %d", &m, &n) != 2)
#else
	if( scanf( "%d/%d", &m, &n) != 2 && scanf( "%d %d", &m, &n) != 2)
#endif
		exit(-1);

	/* Send EvaluatePacket[ Quotient[ m, n]] */
	WSPutFunction( lp, "EvaluatePacket", 1L);
		WSPutFunction( lp, "Quotient", 2L);
			WSPutInteger( lp, m);
			WSPutInteger( lp, n);
	WSEndPacket( lp);
	
	/* skip any packets before the first ReturnPacket */
	while( (pkt = WSNextPacket( lp), pkt) && pkt != RETURNPKT)
		WSNewPacket( lp);
	
	/* inside the ReturnPacket we expect an integer */
	WSGetInteger( lp, &q);
	
	printf( "quotient = %d\n", q);
	
	/* quit Mathematica */
	WSPutFunction( lp, "Exit", 0);

	return 0;
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
