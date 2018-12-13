
#include "wstp.h"
#include <QuEST.h>


QuESTEnv env;

int mytest( int i, int j) {
    
    Qureg q = createQureg(10, env);
    initPlusState(q);
    
    int out = measure(q, 0);
    destroyQureg(q, env);
    return out;
}

qreal anothertest(void) {
    
    Qureg q = createQureg(10, env);
    initPlusState(q);
    
    qreal prob = calcProbOfOutcome(q, 0, 0);
    destroyQureg(q, env);
    return prob;
}


void getQureg(void) {
    
    Qureg q = createQureg(10, env);
    
    
    WSPutFunction(stdlink, "List", 2);
    WSPutInteger(stdlink, 0);
    WSPutString(stdlink, "hello");
    
    
    
    destroyQureg(q, env);
}


void giveQureg(void) {
    
    long int i;
    long int j;
    
    WSGetLongInteger(stdlink, &i);
    WSGetLongInteger(stdlink, &j);
    
    WSPutLongInteger(stdlink, i+j);
}



int main(int argc, char* argv[])
{
    env = createQuESTEnv();
    
	return WSMain(argc, argv);
}
