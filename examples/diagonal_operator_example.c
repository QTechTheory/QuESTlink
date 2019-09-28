#include <stdio.h>
#include "QuEST.h"

int main (int narg, char *varg[]) {
    
    // allocate memory
    int numQubits = 4;
    QuESTEnv env = createQuESTEnv();
    Qureg qureg = createQureg(numQubits, env);
    DiagonalOperator op = createDiagonalOperator(numQubits, env);

    // create input state
    setAmps(qureg, 0, (qreal[16]) {1,2,3,4}, (qreal[16]) {-1,-2,-3,-4}, qureg.numAmpsTotal);
    
    
    printf("state before:\n");
    for (int i=0; i<16; i++) {
        Complex amp = getAmp(qureg, i);
        printf("qureg[%d] = %g + i(%g)\n", i, amp.real, amp.imag);
    }
    
    
    // apply diagonal operator
    for (int i=0; i<8; i++) {
        op.real[i] = .3;
        op.imag[i] = .5;
    }
    syncDiagonalOperator(op);
    applyDiagonalOperator(qureg, op);


    printf("\nstate after:\n");
    for (int i=0; i<16; i++) {
        Complex amp = getAmp(qureg, i);
        printf("qureg[%d] = %g + i(%g)\n", i, amp.real, amp.imag);
    }
    
    
    // free memory
    destroyDiagonalOperator(op);
    destroyQureg(qureg, env); 
    destroyQuESTEnv(env);
    return 0;
}
