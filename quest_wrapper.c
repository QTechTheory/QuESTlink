
#include "wstp.h"
#include <QuEST.h>


QuESTEnv env;

int mytest( int i, int j) {
    
    Qureg q = createQureg(10, env);
    initPlusState(q);
    
    return measure(q, 0);
}

int anothertest(void) {
    
    Qureg q = createQureg(10, env);
    initPlusState(q);
    
    return calcProbOfOutcome(q, 0, 0);
}

int main(int argc, char* argv[])
{
    env = createQuESTEnv();
    
	return WSMain(argc, argv);
}
