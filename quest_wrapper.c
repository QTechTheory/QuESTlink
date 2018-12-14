#include "wstp.h"
#include <QuEST.h>

/**
 * Global instance of QuESTEnv, created when MMA is linked.
 * 
 */
QuESTEnv env;

/**
 * puts a Qureg into MMA, with the structure of
 * {numQubits, isDensityMatrix, realAmps, imagAmps}
 */
void putQuregToMMA(Qureg qureg) {
    
    WSPutFunction(stdlink, "List", 4);
    WSPutInteger(stdlink, qureg.numQubitsRepresented);
    WSPutInteger(stdlink, qureg.isDensityMatrix);
    WSPutReal64List(stdlink, qureg.stateVec.real, qureg.numAmpsTotal);
    WSPutReal64List(stdlink, qureg.stateVec.imag, qureg.numAmpsTotal);
}

/**
 * extracts a Qureg from a List passed from MMA, which
 * must have the structure:
 * {numQubits, isDensityMatrix, realAmps, imagAmps}
 */
Qureg getQuregFromMMA(void) {
    
    // qureg properties
    int numQb, numAmps, isDensityMatrix;
    qreal* realAmps;
    qreal* imagAmps;
    
    // throws error if incorrectly sized arr is passed
    int numElems = 4;
    WSTestHeadWithArgCount(stdlink, "List", &numElems);
    
    // get qureg properties
    WSGetInteger(stdlink, &numQb);
    WSGetInteger(stdlink, &isDensityMatrix);
    WSGetReal64List(stdlink, &realAmps, &numAmps);
    WSGetReal64List(stdlink, &imagAmps, &numAmps);
    
    // validate qurge properties
    // @TODO validate numQb vs numAmps (considering isDensityMatrix)

    // create Qureg
    Qureg qureg;
    if (isDensityMatrix)
        qureg = createDensityQureg(numQb, env);
    else
        qureg = createQureg(numQb, env);
    
    // set wavefunction
    initStateFromAmps(qureg, realAmps, imagAmps);
    
    // free MMA-allocated arrays
    WSReleaseReal64List(stdlink, realAmps, numAmps);
    WSReleaseReal64List(stdlink, imagAmps, numAmps);
    
    return qureg;
}


void giveQureg(void) {
    
    Qureg qureg = getQuregFromMMA();
    qreal ampRe = getRealAmp(qureg, 0);
    
    
    WSPutReal(stdlink, ampRe);
    
    destroyQureg(qureg, env);
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
    
    putQuregToMMA(q);
    
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





int main(int argc, char* argv[])
{
    env = createQuESTEnv();
    
	return WSMain(argc, argv);
}
