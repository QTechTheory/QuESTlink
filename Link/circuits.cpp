/** @file
 * Contains functions for applying circuits to registers.
 *
 * @author Tyson Jones
 */

#include "wstp.h"
#include "QuEST.h"

#include "circuits.hpp"
#include "errors.hpp"
#include "decoders.hpp"
#include "extensions.hpp"
#include "link.hpp"
#include "derivatives.hpp"


/*
 * Max number of target and control qubits which can be specified 
 * for an individual gate 
 */
#define MAX_NUM_TARGS_CTRLS 100

/*
 * Codes for dynamically updating kernel variables, to indicate progress 
 */
#define CIRC_PROGRESS_VAR "QuEST`Private`circuitProgressVar"




int* local_prepareCtrlCache(int* ctrls, int ctrlInd, int numCtrls, int addTarg) {
    static int ctrlCache[MAX_NUM_TARGS_CTRLS]; 
    for (int i=0; i < numCtrls; i++)
        ctrlCache[i] = ctrls[ctrlInd + i];
    if (addTarg != -1)
        ctrlCache[numCtrls] = addTarg;
    return ctrlCache;
}

/* updates the CIRC_PROGRESS_VAR in the front-end with the new passed value 
 * which must lie in [0, 1]. This can be used to indicate progress of a long 
 * evaluation to the user 
 */
void local_updateCircuitProgress(qreal progress) {

    // send new packet to MMA
    WSPutFunction(stdlink, "EvaluatePacket", 1);

    // echo the message
    WSPutFunction(stdlink, "Set", 2);
    WSPutSymbol(stdlink, CIRC_PROGRESS_VAR);
    WSPutReal64(stdlink, progress);

    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
    
    // a new packet is now expected; caller MUST send something else
}





/* @param observablesCache may be NULL
 * @param finalCtrlInd, finalTargInd and finalParamInd are modified to point to 
 *  the final values of ctrlInd, targInd and paramInd, after the #numOps operation 
 *  has been applied. If #numOps isn't smaller than the actual length of the 
 *  circuit arrays, these indices will point out of bounds.
 * @throws QuESTException if a core-QuEST function fails validation (in this case,
 *      exception.thrower will be the name of the throwing core API function),
 *      or if a QuESTlink validation herein fails (exception.thrower will be ""),
 *      or if evaluation is aborted (exception.throw = "Abort")
 */
void local_applyGates(
    Qureg qureg, 
    int numOps, int* opcodes, 
    int* ctrls, int* numCtrlsPerOp, 
    int* targs, int* numTargsPerOp, 
    qreal* params, int* numParamsPerOp,
    qreal* observablesCache,
    int* finalCtrlInd, int* finalTargInd, int* finalParamInd,
    int showProgress
    ) {
        
    int ctrlInd = 0;
    int targInd = 0;
    int paramInd = 0;
    int obvsInd = 0;
    
    // attempt to apply each gate
    for (int opInd=0; opInd < numOps; opInd++) {
                
        // check whether the user has tried to abort
        if (WSMessageReady(stdlink)) {
            int code, arg;
            WSGetMessage(stdlink, &code, &arg);
            if (code == WSTerminateMessage || code == WSInterruptMessage || 
                code == WSAbortMessage     || code == WSImDyingMessage) {
                    
                throw QuESTException("Abort", "Circuit simulation aborted."); // throws
            }
        }
        
        // display progress to the user
        if (showProgress)
            local_updateCircuitProgress(opInd / (qreal) numOps);

        // get gate info
        int op = opcodes[opInd];
        int numCtrls = numCtrlsPerOp[opInd];
        int numTargs = numTargsPerOp[opInd];
        int numParams = numParamsPerOp[opInd];
        
        switch(op) {
            
            case OPCODE_H :
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("Hadamard", numParams, 0); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled Hadamard"); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("Hadamard", numTargs, "1 target"); // throws
                hadamard(qureg, targs[targInd]); // throws
                break;
                
            case OPCODE_S :
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("S gate", numParams, 0); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("S gate", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    sGate(qureg, targs[targInd]); // throws
                else {
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledPhaseShift(qureg, ctrlCache, numCtrls+1, M_PI/2); // throws
                }
                break;
                
            case OPCODE_T : {
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("T gate", numParams, 0); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("T gate", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    tGate(qureg, targs[targInd]); // throws
                else {
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledPhaseShift(qureg, ctrlCache, numCtrls+1, M_PI/4); // throws
                }
            }
                break;
        
            case OPCODE_X : {
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("X", numParams, 0); // throws
                if (numCtrls == 0 && numTargs == 1)
                    pauliX(qureg, targs[targInd]); // throws
                else if (numCtrls == 1 && numTargs == 1)
                    controlledNot(qureg, ctrls[ctrlInd], targs[targInd]); // throws
                else if (numCtrls == 0 && numTargs > 1)
                    multiQubitNot(qureg, &targs[targInd], numTargs); // throws
                else
                    multiControlledMultiQubitNot(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], numTargs); // throws
            }
                break;
                
            case OPCODE_Y :
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("Y", numParams, 0); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("Y", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    pauliY(qureg, targs[targInd]); // throws
                else if (numCtrls == 1)
                    controlledPauliY(qureg, ctrls[ctrlInd], targs[targInd]); // throws
                else
                    throw local_gateUnsupportedExcep("controlled Y"); // throws
                break;
                
            case OPCODE_Z : {
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("Z", numParams, 0); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("Z", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    pauliZ(qureg, targs[targInd]); // throws
                else {
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledPhaseFlip(qureg, ctrlCache, numCtrls+1); // throws
                }
            }
                break;
        
            case OPCODE_Rx :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Rx", numParams, 1); // throws
                if (numCtrls == 0 && numTargs == 1)
                    rotateX(qureg, targs[targInd], params[paramInd]); // throws
                else if (numCtrls == 1 && numTargs == 1)
                    controlledRotateX(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]); // throws
                else {
                    enum pauliOpType paulis[MAX_NUM_TARGS_CTRLS]; 
                    for (int i=0; i<numTargs; i++)
                        paulis[i] = PAULI_X;
                    if (numCtrls == 0)
                        multiRotatePauli(qureg, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
                    else
                        multiControlledMultiRotatePauli(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
                }
                break;
                
            case OPCODE_Ry :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Ry", numParams, 1); // throws
                if (numCtrls == 0 && numTargs == 1)
                    rotateY(qureg, targs[targInd], params[paramInd]); // throws
                else if (numCtrls == 1 && numTargs == 1)
                    controlledRotateY(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]); // throws
                else {
                    enum pauliOpType paulis[MAX_NUM_TARGS_CTRLS]; 
                    for (int i=0; i<numTargs; i++)
                        paulis[i] = PAULI_Y;
                    if (numCtrls == 0)
                        multiRotatePauli(qureg, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
                    else
                        multiControlledMultiRotatePauli(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
                }
                break;
                
            case OPCODE_Rz :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Rz", numParams, 1); // throws
                if (numCtrls == 0 && numTargs == 1)
                    rotateZ(qureg, targs[targInd], params[paramInd]); // throws
                else if (numCtrls == 1 && numTargs == 1)
                    controlledRotateZ(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]); // throws
                else if (numCtrls == 0 && numTargs > 1)
                    multiRotateZ(qureg, &targs[targInd], numTargs, params[paramInd]); // throws
                else
                    multiControlledMultiRotateZ(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], numTargs, params[paramInd]); // throws
                break;
                
            case OPCODE_R: {
                if (numTargs != numParams-1) {
                    throw QuESTException("", 
                        std::string("An internel error in R occured! ") +
                        "The quest_link received an unequal number of Pauli codes " + 
                        "(" + std::to_string(numParams-1) + ") and target qubits " + 
                        "(" + std::to_string(numTargs) + ")!"); // throws
                }
                enum pauliOpType paulis[MAX_NUM_TARGS_CTRLS]; 
                for (int p=0; p < numTargs; p++)
                    paulis[p] = (pauliOpType) ((int) params[paramInd+1+p]);
                if (numCtrls == 0)
                    multiRotatePauli(qureg, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
                else
                    multiControlledMultiRotatePauli(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
            }
                break;
            
            case OPCODE_U : {
                ; // empty post-label statement, courtesy of weird C99 standard.
                // this was added when QuESTlink was already refactored for C++11,
                // however the quirk of being unable to define a variable on the first 
                // line of a switch case in C99 is so strange, I carry on the tradition anyway
                long long int dim = (1 << numTargs);
                if (numParams != 2 * dim*dim)
                    throw QuESTException("", std::to_string(numTargs) + "-qubit U accepts only " + 
                        std::to_string(dim) + "x" +  std::to_string(dim) + " matrices."); // throws
                
                if (numTargs == 1) {
                    ComplexMatrix2 u = local_getMatrix2FromFlatList(&params[paramInd]);
                    if (numCtrls == 0)
                        unitary(qureg, targs[targInd], u); // throws
                    else
                        multiControlledUnitary(qureg, &ctrls[ctrlInd], numCtrls, targs[targInd], u); // throws
                }
                else if (numTargs == 2) {
                    ComplexMatrix4 u = local_getMatrix4FromFlatList(&params[paramInd]);
                    if (numCtrls == 0)
                        twoQubitUnitary(qureg, targs[targInd], targs[targInd+1], u); // throws
                    else
                        multiControlledTwoQubitUnitary(qureg, &ctrls[ctrlInd], numCtrls, targs[targInd], targs[targInd+1], u); // throws
                } 
                else {
                    // this is wastefully(?) allocating and deallocating memory on the fly!
                    ComplexMatrixN u = createComplexMatrixN(numTargs);
                    local_setMatrixNFromFlatList(&params[paramInd], u, numTargs);
                    if (numCtrls == 0)
                        multiQubitUnitary(qureg, &targs[targInd], numTargs, u); // throws
                    else
                        multiControlledMultiQubitUnitary(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], numTargs, u); // throws
                    // memory leak if above throws :^)
                    destroyComplexMatrixN(u);
                }
            }
                break;
                
            case OPCODE_Matr : {
                ;
                long long int dim = (1 << numTargs);
                if (numParams != 2 * dim*dim)
                    throw QuESTException("", std::to_string(numTargs) + "-qubit Matr accepts only " + 
                        std::to_string(dim) + "x" +  std::to_string(dim) + " matrices."); // throws
                
                // this is wastefully(?) allocating and deallocating memory on the fly!
                ComplexMatrixN m = createComplexMatrixN(numTargs);
                local_setMatrixNFromFlatList(&params[paramInd], m, numTargs);
                if (numCtrls == 0)
                    applyGateMatrixN(qureg, &targs[targInd], numTargs, m); // throws
                else
                    applyMultiControlledGateMatrixN(qureg, &ctrls[ctrlInd], numCtrls, &targs[targInd], numTargs, m); // throws
                // memory leak if above throws :^)
                destroyComplexMatrixN(m);
                
            }
                break;
                
            case OPCODE_Deph :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Dephasing", numParams, 1); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled dephasing"); // throws
                if (numTargs != 1 && numTargs != 2)
                    throw local_wrongNumGateTargsExcep("Dephasing", numTargs, "1 or 2 targets"); // throws
                if (params[paramInd] == 0)
                    break; // allows zero-prob decoherence to act on state-vectors
                if (numTargs == 1)
                    mixDephasing(qureg, targs[targInd], params[paramInd]); // throws
                if (numTargs == 2)
                    mixTwoQubitDephasing(qureg, targs[targInd], targs[targInd+1], params[paramInd]); // throws
                break;
                
            case OPCODE_Depol :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Depolarising", numParams, 1); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled depolarising"); // throws
                if (numTargs != 1 && numTargs != 2)
                    throw local_wrongNumGateTargsExcep("Depolarising", numTargs, "1 or 2 targets"); // throws
                if (params[paramInd] == 0)
                    break; // allows zero-prob decoherence to act on state-vectors
                if (numTargs == 1)
                    mixDepolarising(qureg, targs[targInd], params[paramInd]); // throws
                if (numTargs == 2)
                    mixTwoQubitDepolarising(qureg, targs[targInd], targs[targInd+1], params[paramInd]); // throws
                break;
                
            case OPCODE_Damp :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Damping", numParams, 1); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled damping"); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("Damping", numTargs, "1 target"); // throws
                if (params[paramInd] == 0)
                    break; // allows zero-prob decoherence to act on state-vectors
                mixDamping(qureg, targs[targInd], params[paramInd]); // throws
                break;
                
            case OPCODE_SWAP: {
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("SWAP", numParams, 0); // throws
                if (numTargs != 2)
                    throw local_wrongNumGateTargsExcep("Depolarising", numTargs, "2 targets"); // throws
                if (numCtrls == 0) {
                    swapGate(qureg, targs[targInd],   targs[targInd+1]); // throws
                } else {    
                    // core-QuEST doesn't yet support multiControlledSwapGate, 
                    // so we construct SWAP from 3 CNOT's, and add additional controls
                    ComplexMatrix2 u;
                    u.real[0][0] = 0; u.real[0][1] = 1; // verbose for old MSVC 
                    u.real[1][0] = 1; u.real[1][1] = 0;
                    u.imag[0][0] = 0; u.imag[0][1] = 0;
                    u.imag[1][0] = 0; u.imag[1][1] = 0;
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledUnitary(qureg, ctrlCache, numCtrls+1, targs[targInd+1], u); // throws
                    ctrlCache[numCtrls] = targs[targInd+1];
                    multiControlledUnitary(qureg, ctrlCache, numCtrls+1, targs[targInd], u);
                    ctrlCache[numCtrls] = targs[targInd];
                    multiControlledUnitary(qureg, ctrlCache, numCtrls+1, targs[targInd+1], u);
                }
            }
                break;
                
            case OPCODE_M: {
                if (numParams != 0)
                    throw local_wrongNumGateParamsExcep("M", numParams, 0); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled measurement"); // throws
                for (int q=0; q < numTargs; q++) {
                    int outcomeVal = measure(qureg, targs[targInd+q]); // throws // throws
                    if (observablesCache != NULL)
                        observablesCache[obvsInd++] = outcomeVal;
                }
            }
                break;
            
            case OPCODE_P:
                if (numParams != 1 && numParams != numTargs)
                    throw QuESTException("", 
                        std::string("P[outcomes] specified a different number of binary outcomes ") + 
                        "(" + std::to_string(numParams) + ") than target qubits  " +
                        "(" + std::to_string(numTargs) + ")!"); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled projector"); // throws
                if (numParams > 1) {
                    qreal prob = 1;
                    for (int q=0; q < numParams; q++)
                        prob *= collapseToOutcome(qureg, targs[targInd+q], (int) params[paramInd+q]); // throws
                    if (observablesCache != NULL)
                        observablesCache[obvsInd++] = prob;
                }
                else {
                    // check value isn't impossibly high
                    if (params[paramInd] >= (1LL << numTargs))
                        throw QuESTException("",
                            "P[ " + std::to_string((int) params[paramInd]) + "] was applied to " +
                            std::to_string(numTargs) + " qubits and exceeds their maximum represented " +
                            "value of " + std::to_string(1LL << numTargs) + "."); // throws
                    // work out each bit outcome and apply; right most (least significant) bit acts on right-most target
                    qreal prob = 1;
                    for (int q=0; q < numTargs; q++)
                        prob *= collapseToOutcome(qureg, targs[targInd+numTargs-q-1], (((int) params[paramInd]) >> q) & 1); // throws
                    if (observablesCache != NULL)
                        observablesCache[obvsInd++] = prob;
                }
                break;
                
            case OPCODE_Kraus: {
                ; // empty post-label statement, courtesy of weird C99 standard
                int numKrausOps = (int) params[paramInd];
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled Kraus map"); // throws
                if (numTargs != 1 && numTargs != 2)
                    throw local_wrongNumGateTargsExcep("Kraus map", numTargs, "1 or 2 targets"); // throws
                if ((numKrausOps < 1) ||
                    (numTargs == 1 && numKrausOps > 4 ) ||
                    (numTargs == 2 && numKrausOps > 16))
                    throw QuESTException("", 
                        std::to_string(numKrausOps) + " operators were passed to " +
                        std::to_string(numTargs) +  "-qubit Kraus[ops], which accepts only >0 and <=" + 
                        std::to_string((numTargs==1)? 4:16) + " operators!"); // throws
                if (numTargs == 1 && (numParams-1) != 2*2*2*numKrausOps)
                    throw QuESTException("", "one-qubit Kraus expects 2-by-2 matrices!"); // throws
                if (numTargs == 2 && (numParams-1) != 4*4*2*numKrausOps)
                    throw QuESTException("", "two-qubit Kraus expects 4-by-4 matrices!"); // throws

                if (numTargs == 1) {
                    ComplexMatrix2 krausOps[4];
                    int opElemInd = 1 + paramInd;
                    for (int n=0; n < numKrausOps; n++)
                        krausOps[n] = local_getMatrix2FromFlatList(&params[opElemInd + 2*2*2*n]);
                    mixKrausMap(qureg, targs[targInd], krausOps, numKrausOps); // throws
                } 
                else if (numTargs == 2) {
                    ComplexMatrix4 krausOps[16];
                    int opElemInd = 1 + paramInd;
                    for (int n=0; n < numKrausOps; n++)
                        krausOps[n] = local_getMatrix4FromFlatList(&params[opElemInd + 2*4*4*n]);
                    mixTwoQubitKrausMap(qureg, targs[targInd], targs[targInd+1], krausOps, numKrausOps); // throws
                }
            }
                break;
                
            case OPCODE_KrausNonTP: {
                ; // empty post-label statement, courtesy of weird C99 standard
                int numKrausOps = (int) params[paramInd];
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled non-trace-preserving Kraus map"); // throws
                if (numTargs != 1 && numTargs != 2)
                    throw local_wrongNumGateTargsExcep("non-trace-preserving Kraus map", numTargs, "1 or 2 targets"); // throws
                if ((numKrausOps < 1) ||
                    (numTargs == 1 && numKrausOps > 4 ) ||
                    (numTargs == 2 && numKrausOps > 16))
                    throw QuESTException("", 
                        std::to_string(numKrausOps) + " operators were passed to " +
                        std::to_string(numTargs) +  "-qubit KrausNonTP[ops], which accepts only >0 and <=" + 
                        std::to_string((numTargs==1)? 4:16) + " operators!"); // throws
                if (numTargs == 1 && (numParams-1) != 2*2*2*numKrausOps)
                    throw QuESTException("", "one-qubit non-trace-preserving Kraus expects 2-by-2 matrices!"); // throws
                if (numTargs == 2 && (numParams-1) != 4*4*2*numKrausOps)
                    throw QuESTException("", "two-qubit non-trace-preserving Kraus expects 4-by-4 matrices!"); // throws

                if (numTargs == 1) {
                    ComplexMatrix2 krausOps[4];
                    int opElemInd = 1 + paramInd;
                    for (int n=0; n < numKrausOps; n++)
                        krausOps[n] = local_getMatrix2FromFlatList(&params[opElemInd + 2*2*2*n]);
                    mixNonTPKrausMap(qureg, targs[targInd], krausOps, numKrausOps); // throws
                } 
                else if (numTargs == 2) {
                    ComplexMatrix4 krausOps[16];
                    int opElemInd = 1 + paramInd;
                    for (int n=0; n < numKrausOps; n++)
                        krausOps[n] = local_getMatrix4FromFlatList(&params[opElemInd + 2*4*4*n]);
                    mixNonTPTwoQubitKrausMap(qureg, targs[targInd], targs[targInd+1], krausOps, numKrausOps); // throws
                }
            }
                break;
                
            case OPCODE_G : {
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Global phase", numParams, 1); // throws
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled global phase"); // throws
                if (numTargs != 0)
                    throw local_wrongNumGateTargsExcep("Global phase", numTargs, "0 targets"); // throws
                if (params[paramInd] == 0)
                    break;
                // phase does not change density matrices
                if (!qureg.isDensityMatrix) {
                     // create factor exp(i param)
                    Complex zero; zero.real=0; zero.imag=0;
                    Complex fac; fac.real=cos(params[paramInd]); fac.imag=sin(params[paramInd]);
                    setWeightedQureg(zero, qureg, zero, qureg, fac, qureg); // throws
                }
            }
                break;
                
            case OPCODE_Id :
                // any numCtrls, numParams and numTargs is valid; all do nothing!
                break;
                
            case OPCODE_Ph : {
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Ph", numParams, 1); // throws
                int numQubits = numCtrls + numTargs;
                if (numQubits < 1)
                    throw local_wrongNumGateTargsExcep("Ph", numQubits, "at least 1 qubit (between control and target qubits)"); // throws
                if (params[paramInd] == 0)
                    break;
                    
                // unpack all controls and targets (since symmetric)
                // (ctrlCache has length [MAX_NUM_TARGS_CTRLS], so it can fit all targs)
                int* qubitCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, -1);
                for (int i=0; i<numTargs; i++)
                    qubitCache[numCtrls+i] = targs[targInd+i];
                
                // attempt optimisations first
                if (numQubits == 1)
                    phaseShift(qureg, qubitCache[0], params[paramInd]);
                else if (numQubits == 2)
                    controlledPhaseShift(qureg, qubitCache[0], qubitCache[1], params[paramInd]);
                else
                    multiControlledPhaseShift(qureg, qubitCache, numQubits, params[paramInd]);
            }
                break;
                
            default:            
                throw QuESTException("", "circuit contained an unknown gate."); // throws
        }
        
        // update progress through gate list
        ctrlInd += numCtrls;
        targInd += numTargs;
        paramInd += numParams;
    }
    
    // update final pointers
    *finalCtrlInd = ctrlInd;
    *finalTargInd = targInd;
    *finalParamInd = paramInd;
}


/* Applies a given circuit to the identified qureg.
 * The circuit is expressed as lists of opcodes (identifying gates),
 * the total flat sequence control qubits, a list denoting how many of
 * the control qubits apply to each operation, their target qubits (flat list),
 * a list denotating how many targets each operation has, their parameters 
 * (flat list) and a list denoting how many params each operation has.
 * Returns a list of measurement outcome of any performed measurements in the circuit.
 * The original qureg of the state is restored when this function
 * is aborted by the calling MMA (sends Abort[] to MMA), or aborted due to encountering
 * an invalid gate or a QuEST-core validation error (sends $Failed to MMA). 
 */
void internal_applyCircuit(int id, int storeBackup, int showProgress) {
    
    // get arguments from MMA link; these must be later freed!
    int numOps;
    int *opcodes, *ctrls, *numCtrlsPerOp, 
        *targs, *numTargsPerOp, *numParamsPerOp;
    qreal* params;
    int totalNumCtrls, totalNumTargs, totalNumParams; // these fields are only needed by clean-up
    local_loadCircuitFromMMA(
        &numOps, &opcodes, &ctrls, &numCtrlsPerOp, 
        &targs, &numTargsPerOp, &params, &numParamsPerOp,
        &totalNumCtrls, &totalNumTargs, &totalNumParams);
    
    // ensure qureg exists, else clean-up and exit
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        
    } catch (QuESTException& err) {
    
        local_freeCircuit(
            opcodes, ctrls, numCtrlsPerOp, targs, 
            numTargsPerOp, params, numParamsPerOp,
            numOps, totalNumCtrls, totalNumTargs, totalNumParams);
            
        local_sendErrorAndFail("ApplyCircuit", err.message);
        return;
    }
    
    Qureg qureg = quregs[id];
    Qureg backup;
    if (storeBackup)
        backup = createCloneQureg(qureg, env); // must clean-up
    
    // count the total number of measurements performed in a circuit
    int totalNumObvsGates = 0;
    int totalNumObvs = 0;
    for (int opInd=0; opInd < numOps; opInd++) {
        if (opcodes[opInd] == OPCODE_M || opcodes[opInd] == OPCODE_P) {
            totalNumObvsGates++;
            totalNumObvs += (opcodes[opInd] == OPCODE_P)? 1 : numTargsPerOp[opInd];
        }
    }
        
    // prepare records of measurement outcomes
    qreal* observablesCache = (qreal*) malloc(totalNumObvs * sizeof(qreal)); // must clean-up
    int obsInd = 0;

    // attempt to apply the circuit
    try {
        // these fields are ignored
        int finalCtrlInd, finalTargInd, finalParamInd;
        
        local_applyGates(
            qureg, numOps, opcodes, ctrls, numCtrlsPerOp, 
            targs, numTargsPerOp, params, numParamsPerOp,
            observablesCache,
            &finalCtrlInd, &finalTargInd, &finalParamInd,
            showProgress); // throws
            
        // return lists of measurement outcomes
        obsInd = 0;
        WSPutFunction(stdlink, "List", totalNumObvsGates);
        for (int opInd=0; opInd < numOps; opInd++) {
            if (opcodes[opInd] == OPCODE_M) {
                // M_0,1,2 -> {0, 1, 1}, M_0 -> {1}
                WSPutFunction(stdlink, "List", numTargsPerOp[opInd]);
                for (int i=0; i < numTargsPerOp[opInd]; i++)
                    WSPutInteger(stdlink, (int) observablesCache[obsInd++]);
            }
            if (opcodes[opInd] == OPCODE_P) {
                // P_0 -> .1
                WSPutReal64(stdlink, observablesCache[obsInd++]);
            }
        }
        
        // clean-up
        if (storeBackup)
            destroyQureg(backup, env);
        free(observablesCache);
        local_freeCircuit(
            opcodes, ctrls, numCtrlsPerOp, targs, 
            numTargsPerOp, params, numParamsPerOp,
            numOps, totalNumCtrls, totalNumTargs, totalNumParams);
    
    } catch (QuESTException& err) {
        
        // restore backup (if made)
        if (storeBackup) {
            cloneQureg(qureg, backup);
            destroyQureg(backup, env);
        }
                
        // all objs need cleaning
        free(observablesCache);
        local_freeCircuit(
            opcodes, ctrls, numCtrlsPerOp, targs, 
            numTargsPerOp, params, numParamsPerOp,
            numOps, totalNumCtrls, totalNumTargs, totalNumParams);
            
        // report error, depending on type
        std::string backupNotice;
        if (storeBackup)
             backupNotice = " The qureg (id " + std::to_string(id) + 
                ") has been restored to its prior state.";
        else
            backupNotice = " Since no backup was stored, the qureg (id " + std::to_string(id) + 
                ") is now in an unknown state, and should be reinitialised.";
        
        if (err.thrower == "")
            local_sendErrorAndFail("ApplyCircuit", err.message + backupNotice);
        else if (err.thrower == "Abort")
            local_sendErrorAndAbort("ApplyCircuit", err.message + backupNotice);
        else 
            local_sendErrorAndFail("ApplyCircuit", 
                "Error in " + err.thrower + ": " + err.message + backupNotice);
    }
}





