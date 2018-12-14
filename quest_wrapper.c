
#include "wstp.h"
#include <QuEST.h>


QuESTEnv env;


/**
 * puts Rule[field, value] in current MMA expression
 */
void putRuleInt(char* field, int value) {
    
    WSPutFunction(stdlink, "Rule", 2);
    WSPutString(stdlink, field);
    WSPutInteger(stdlink, value);
}

/**
 * puts Rule[field, arr] in current MMA expression
 */
void putRuleArr(char* field, qreal* arr, int len) {
    
    WSPutFunction(stdlink, "Rule", 2);
    WSPutString(stdlink, field);
    WSPutRealList(stdlink, arr, len);
}

/**
 * puts an Association representing a Qureg in current MMA expression
 */
void putQureg(Qureg q) {

    WSPutFunction(stdlink, "Association", 6);
    putRuleInt("numQubitsRepresented", q.numQubitsRepresented);
    putRuleInt("numQubitsInStateVec", q.numQubitsInStateVec);
    putRuleInt("isDensityMatrix", q.isDensityMatrix);
    putRuleInt("numAmpsTotal", q.numAmpsTotal);
    putRuleArr("stateVecReal", q.stateVec.real, q.numAmpsTotal);
    putRuleArr("stateVecImag", q.stateVec.imag, q.numAmpsTotal);
}

/**
 * extracts a Qureg from an Assocation passed from MMA
 */
void getQureg(void) {
    
    // throws error if not a correctly sized Assocation was passed
    int numArgs = 6;
    WSTestHeadWithArgCount(stdlink, "Association", &numArgs);

    

}



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


void returnQureg(void) {
    
    Qureg q = createQureg(4, env);
    initPlusState(q);
    controlledNot(q, 0, 1);
    controlledPhaseShift(q, 1, 2, .3);
    
    putQureg(q);
    
    destroyQureg(q, env);
    
    // WSPutRealList
    // WSPutDoubleArray
    
    
    // for density matrices:
    
    // mat = (double *) calloc(n * m, sizeof(double))
    // mat[i + j*n] = ...
    // long dimensions[2]
    // dimensions[0] = n; dimensions[1] = m;
    // int depth = 2;
    // char *heads[2] = {"List", "List"}
    // MLPutDoubleArray(stdlink, mat, dimensions, heads, depth);
    // free(mat)
}


void giveQureg(void) {
    
    getQureg();
    
    

}



int main(int argc, char* argv[])
{
    env = createQuESTEnv();
    
	return WSMain(argc, argv);
}
