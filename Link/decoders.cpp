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



/*
 * parsers
 */
 
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



/*
 * Circuit methods 
 */
 
void Circuit::loadFromMMA() {
    
    // prepare pointers (left will persist, middle won't)
    int* opcodes;
    int* ctrls;     int* numCtrlsPerOp;
    int* targs;     int* numTargsPerOp;
    qreal* params;  int* numParamsPerOp;
    
    // load arrays from MMA (order fixed by QuESTlink.m)
    WSGetInteger32List(stdlink, &opcodes, &numGates);
    WSGetInteger32List(stdlink, &ctrls, &totalNumCtrls);
    WSGetInteger32List(stdlink, &numCtrlsPerOp, &numGates);
    WSGetInteger32List(stdlink, &targs, &totalNumTargs);
    WSGetInteger32List(stdlink, &numTargsPerOp, &numGates);
    WSGetReal64List(stdlink, &params, &totalNumParams);
    WSGetInteger32List(stdlink, &numParamsPerOp, &numGates);
    
    // allocate gates array attribute (creates all Gate instances)
    gates = new Gate[numGates];
    
    int ctrlInd = 0;
    int targInd = 0;
    int paramInd = 0;
    
    for (int opInd=0; opInd<numGates; opInd++) {
        
        int numCtrls  = numCtrlsPerOp[opInd];
        int numTargs  = numTargsPerOp[opInd];
        int numParams = numParamsPerOp[opInd];
        
        // initialise each gate, persisting ctrls, targs, params
        gates[opInd].init(
            opcodes[opInd], 
            &ctrls[ctrlInd],   numCtrls,
            &targs[targInd],   numTargs,
            &params[paramInd], numParams);
            
        ctrlInd  += numCtrls;
        targInd  += numTargs;
        paramInd += numParams;
    }
    
    // clean-up the arrays which are no longer needed
    WSReleaseInteger32List(stdlink, opcodes, numGates);
    WSReleaseInteger32List(stdlink, numCtrlsPerOp,  numGates);
    WSReleaseInteger32List(stdlink, numTargsPerOp,  numGates);
    WSReleaseInteger32List(stdlink, numParamsPerOp, numGates);
}

void Circuit::freeMMA() {
    
    // first gate has circuit-wide gate array pointers
    Gate gate = gates[0];
    WSReleaseInteger32List(stdlink, gate.getCtrlsAddr(),  totalNumCtrls);
    WSReleaseInteger32List(stdlink, gate.getTargsAddr(),  totalNumTargs);
    WSReleaseReal64List(   stdlink, gate.getParamsAddr(), totalNumParams);
}

void Circuit::sendOutputsToMMA(qreal* outputs) {
    
    // let each output-gate correspond to an element or a sublist of the MMA array
    WSPutFunction(stdlink, "List", getNumGatesWithOutputs());
    
    int outInd = 0;
    
    for (int gateInd=0; gateInd<numGates; gateInd++) {
        
        Gate gate = gates[gateInd];
        int numOuts = gate.getNumOutputs();
        
        // P gate outputs are top-level qreals (probabilities of projector)
        if (gate.getOpcode() == OPCODE_P)
            WSPutReal64(stdlink, outputs[outInd++]);
        
        // M gate outputs are integer sublists (measurement outcomes grouped by targets)
        if (gate.getOpcode() == OPCODE_M) {
            WSPutFunction(stdlink, "List", numOuts);
            for (int i=0; i<numOuts; i++)
                WSPutInteger(stdlink, (int) outputs[outInd++]);
        }
    }
}



/*
 * Hamiltonian loading
 */

void local_loadEncodedPauliSumFromMMA(int* numPaulis, int* numTerms, qreal** termCoeffs, int** allPauliCodes, int** allPauliTargets, int** numPaulisPerTerm) {
    
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
