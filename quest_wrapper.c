
#include "wstp.h"
#include <QuEST.h>


QuESTEnv env;


extern int WSMain(int, char **);



extern int mytest( int i, int j);

int mytest( int i, int j)
{
    
    Qureg q = createQureg(10, env);
    initPlusState(q);
    
    return measure(q, 0);
    
    
	//return i+j;
}

int main(int argc, char* argv[])
{
    env = createQuESTEnv();
    
	return WSMain(argc, argv);
}
