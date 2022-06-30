/* @file
 * Translates between Mathematica encodings of structures like matrices and circuits, 
 * and their C++ and/or representations.
 *
 * @author Tyson Jones
 */

#include "QuEST.h"
#include "wstp.h"

#include "errors.hpp"
#include "circuits.hpp"

#include <string>


ComplexMatrix2 local_getMatrix2FromFlatList(qreal* list) {
    int dim = 2;
    ComplexMatrix2 m;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            m.real[r][c] = list[2*(dim*r+c)];
            m.imag[r][c] = list[2*(dim*r+c)+1];
        }
    return m;
}

ComplexMatrix4 local_getMatrix4FromFlatList(qreal* list) {
    int dim = 4;
    ComplexMatrix4 m;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            m.real[r][c] = list[2*(dim*r+c)];
            m.imag[r][c] = list[2*(dim*r+c)+1];
        }
    return m;
}

void local_setMatrixNFromFlatList(qreal* list, ComplexMatrixN m, int numQubits) {
    long long int dim = (1LL << numQubits);
    for (long long int r=0; r<dim; r++)
        for (long long int c=0; c<dim; c++) {
            m.real[r][c] = list[2*(dim*r+c)];
            m.imag[r][c] = list[2*(dim*r+c)+1];
        }
}




/* load a circuit specification from Mathematica. All these lists must be 
 * later freed
 */
void local_loadCircuitFromMMA(
    int* numOps, int** opcodes, int** ctrls, int** numCtrlsPerOp, 
    int** targs, int** numTargsPerOp, qreal** params, int** numParamsPerOp,
    int* totalNumCtrls, int* totalNumTargs, int* totalNumParams
) {
    WSGetInteger32List(stdlink, opcodes, numOps);
    
    WSGetInteger32List(stdlink, ctrls, totalNumCtrls);
    WSGetInteger32List(stdlink, numCtrlsPerOp, numOps);

    WSGetInteger32List(stdlink, targs, totalNumTargs);
    WSGetInteger32List(stdlink, numTargsPerOp, numOps);

    WSGetReal64List(stdlink, params, totalNumParams);
    WSGetInteger32List(stdlink, numParamsPerOp, numOps);
}

void local_freeCircuit(
    int* opcodes, int* ctrls, int* numCtrlsPerOp, int* targs, 
    int* numTargsPerOp, qreal* params, int* numParamsPerOp,
    int numOps, int totalNumCtrls, int totalNumTargs, int totalNumParams
) {
    WSReleaseInteger32List(stdlink, opcodes, numOps);
    
    WSReleaseInteger32List(stdlink, ctrls, totalNumCtrls);
    WSReleaseInteger32List(stdlink, numCtrlsPerOp, numOps);
    
    WSReleaseInteger32List(stdlink, targs, totalNumTargs);
    WSReleaseInteger32List(stdlink, numTargsPerOp, numOps);
    
    WSReleaseReal64List(stdlink, params, totalNumParams);
    WSReleaseInteger32List(stdlink, numParamsPerOp, numOps);
}




void local_loadEncodedPauliSumFromMMA(int* numPaulis, int* numTerms, qreal** termCoeffs, int** allPauliCodes, int** allPauliTargets, int** numPaulisPerTerm) {
    // get arguments from MMA link
    WSGetReal64List(stdlink, termCoeffs, numTerms);
    WSGetInteger32List(stdlink, allPauliCodes, numPaulis);    
    WSGetInteger32List(stdlink, allPauliTargets, numPaulis);
    WSGetInteger32List(stdlink, numPaulisPerTerm, numTerms);
}

/* can throw exception if pauli targets are invalid indices.
 * in that event, arrPaulis will be freed internally, and set to NULL (if init'd by caller)
 */
pauliOpType* local_decodePauliSum(int numQb, int numTerms, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm) {
    
    // convert {allPauliCodes}, {allPauliTargets}, {numPaulisPerTerm}, and
    // qureg.numQubitsRepresented into {pauli-code-for-every-qubit}
    int arrLen = numTerms * numQb;
    pauliOpType* arrPaulis = (pauliOpType*) malloc(arrLen * sizeof *arrPaulis);
    for (int i=0; i < arrLen; i++)
        arrPaulis[i] = PAULI_I;
    
    int allPaulisInd = 0;
    for (int t=0;  t < numTerms; t++) {
        for (int j=0; j < numPaulisPerTerm[t]; j++) {
            int currTarget = allPauliTargets[allPaulisInd];
            
            // clean-up and error on invalid target qubit
            if (currTarget >= numQb) {
                free(arrPaulis);
                arrPaulis = NULL;
                throw QuESTException("",
                    "Invalid target index (" + std::to_string(currTarget) + 
                    ") of Pauli operator in Pauli sum of " + std::to_string(numQb) + " qubits.");
                }
            
            int arrInd = t*numQb + currTarget;
            int code = allPauliCodes[allPaulisInd++];
            pauliOpType op = (code == OPCODE_Id)? PAULI_I : (pauliOpType) code;
            arrPaulis[arrInd] = op;
        }
    }
    
    // must be freed by caller
    return arrPaulis;
}

void local_freePauliSum(int numPaulis, int numTerms, qreal* termCoeffs, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm, pauliOpType* arrPaulis) {
    WSReleaseReal64List(stdlink, termCoeffs, numTerms);
    WSReleaseInteger32List(stdlink, allPauliCodes, numPaulis);
    WSReleaseInteger32List(stdlink, allPauliTargets, numPaulis);
    WSReleaseInteger32List(stdlink, numPaulisPerTerm, numTerms);
    
    // may be none if validation triggers clean-up before arrPaulis gets created
    if (arrPaulis != NULL)
        free(arrPaulis);
}
