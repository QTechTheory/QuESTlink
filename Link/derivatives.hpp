
#ifndef DERIVATIVES_H
#define DERIVATIVES_H


void local_getDerivativeQuregs(
    // variable (to be differentiated) info
    int* quregIds, int* varOpInds, int numVars, 
    // circuit info
    int numOps, int* opcodes, 
    int* ctrls, int* numCtrlsPerOp, 
    int* targs, int* numTargsPerOp, 
    qreal* params, int* numParamsPerOp,
    // derivative matrices of general unitary gates in circuit
    qreal* unitaryDerivs);
    
void internal_calcQuregDerivs(int initStateId);


#endif // DERIVATIVES_H