/** @file 
 * Contains functions for processing analytic derivatives of circuits and channels.
 */

#include "wstp.h"

#include "QuEST.h"
#include "QuEST_internal.h"
#include "QuEST_cpu_internal.h"
#include "QuEST_validation.h"

#include "errors.hpp"
#include "decoders.hpp"
#include "extensions.hpp"
#include "circuits.hpp"
#include "link.hpp"



/* @precondition quregs must be prior initialised and cloned to the initial state of the circuit.
 * @throws QuESTException if a core QuEST validation fails (in this case,
 *      exception.thrower will be the name of the throwing core API function), 
 *      or we encounter an unsupported gate (exception.thrower = ""), 
 *      or if evaluation is aborted (exception.throw = "Abort").
 *      If an exception is thrown, some quregs may have been modified and are not restored.
 */
void local_getDerivativeQuregs(
    // variable (to be differentiated) info
    int* quregIds, int* varOpInds, int numVars, 
    // circuit info
    int numOps, int* opcodes, 
    int* ctrls, int* numCtrlsPerOp, 
    int* targs, int* numTargsPerOp, 
    qreal* params, int* numParamsPerOp,
    // derivative matrices of general unitary gates in circuit
    qreal* unitaryDerivs)
{        
    // index of the first (real) element of the next unitary derivative
    int unitaryDerivInd = 0;

    // don't record observables
    qreal* observables = NULL;
    
    // don't dynamically update frontend with progress
    int dontShowProgress = 0;
    
    // compute each derivative one-by-one
    for (int v=0; v<numVars; v++) {
        
        // check whether the user has tried to abort
        if (WSMessageReady(stdlink)) {
            int code, arg;
            WSGetMessage(stdlink, &code, &arg);
            if (code == WSTerminateMessage || code == WSInterruptMessage || 
                code == WSAbortMessage     || code == WSImDyingMessage) {
                    
                throw QuESTException("Abort", "Circuit simulation aborted."); // throws
            }
        }
        
        local_throwExcepIfQuregNotCreated(quregIds[v]); // throws
        Qureg qureg = quregs[quregIds[v]];
        int varOp = varOpInds[v];
        
        // indices AFTER last gate applied by circuit
        int finalCtrlInd, finalTargInd, finalParamInd;
        
        // apply only gates up to and including the to-be-differentiated gate,
        // unless that gate is the general unitary
        int op = opcodes[varOp];
        int diffGateWasApplied = (op != OPCODE_U);
        local_applyGates(
            qureg, (diffGateWasApplied)? varOp+1 : varOp, opcodes, 
            ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp,
            observables, &finalCtrlInd, &finalTargInd, &finalParamInd, 
            dontShowProgress); // throws 

        // details of (possibly already applied) to-be-differentiated gate
        int numCtrls = numCtrlsPerOp[varOp];
        int numTargs = numTargsPerOp[varOp];
        int numParams = numParamsPerOp[varOp];

        // wind back inds to point to the to-be-differentiated gate 
        if (diffGateWasApplied) {
            finalCtrlInd -= numCtrls;
            finalTargInd -= numTargs;
            finalParamInd -= numParams;
        }
        
        // choices of re-normalisation (verbose for MSVC :( )
        Complex negHalfI; negHalfI.real=0; negHalfI.imag=-0.5;
        Complex posI; posI.real=0; posI.imag=1;
        Complex zero; zero.real=0; zero.imag=0;
        Complex one; one.real=1; one.imag=0;
        
        // disregard control qubits and apply gate Paulis incurred by differentiation 
        Complex normFac;
        switch(op) {
            case OPCODE_Rx:
                for (int t=0; t < numTargs; t++) // multi-target X may be possible later 
                    pauliX(qureg, targs[t+finalTargInd]);  // throws
                normFac = negHalfI;
                break;
            case OPCODE_Ry:
                for (int t=0; t < numTargs; t++) // multi-target Y may be possible later 
                    pauliY(qureg, targs[t+finalTargInd]);  // throws
                normFac = negHalfI;
                break;
            case OPCODE_Rz:
                for (int t=0; t < numTargs; t++)
                    pauliZ(qureg, targs[t+finalTargInd]);  // throws
                normFac = negHalfI;
                break;
            case OPCODE_R:
                for (int t=0; t < numTargs; t++) {
                    int pauliCode = (int) params[t+finalParamInd+1];
                    if (pauliCode == 1) pauliX(qureg, targs[t+finalTargInd]);  // throws
                    if (pauliCode == 2) pauliY(qureg, targs[t+finalTargInd]);  // throws
                    if (pauliCode == 3) pauliZ(qureg, targs[t+finalTargInd]);  // throws
                }
                normFac = negHalfI;
                break;
            case OPCODE_G:
                ; // no additional gate introduced by derivative
                normFac = posI;
                break;
            case OPCODE_U:
                if (numTargs == 1) {
                    ComplexMatrix2 u2 = local_getMatrix2FromFlatList(&unitaryDerivs[unitaryDerivInd]);
                    unitaryDerivInd += 2*2*2;
                    applyMatrix2(qureg, targs[finalTargInd], u2); // throws
                } else if (numTargs == 2) {
                    ComplexMatrix4 u4 = local_getMatrix4FromFlatList(&unitaryDerivs[unitaryDerivInd]);
                    unitaryDerivInd += 2*4*4;
                    applyMatrix4(qureg, targs[finalTargInd], targs[finalTargInd+1], u4); // throws
                }
                else {
                    // TODO: create a non-dynamic ComplexMatrixN instance 
                    // general matrix N unpacking here; can do static
                    throw QuESTException("", "multi-qubit U derivative is not yet supported."); // throws
                }
                normFac = one;
                break;
            default:            
                throw QuESTException("", "Only Rx, Ry, Rz, R, U and their controlled gates may be differentiated."); // throws
        }
        
        // differentiate control qubits by forcing them to 1, without renormalising
        for (int c=0; c<numCtrls; c++)
            applyProjector(qureg, ctrls[finalCtrlInd], 1); // throws
        
        // adjust normalisation
        setWeightedQureg(zero, qureg, zero, qureg, normFac, qureg); // cannot throw
        
        // wind forward inds to point to the next gate 
        finalCtrlInd += numCtrls;
        finalTargInd += numTargs;
        finalParamInd += numParams;

        // apply the remainder of the circuit
        local_applyGates(
            qureg, numOps-(varOp+1), &opcodes[varOp+1], 
            &ctrls[  finalCtrlInd], &numCtrlsPerOp[varOp+1], 
            &targs[  finalTargInd], &numTargsPerOp[varOp+1], 
            &params[finalParamInd], &numParamsPerOp[varOp+1],
            observables, &finalCtrlInd, &finalTargInd, &finalParamInd, 
            dontShowProgress); // throws
    }
}

void internal_calcQuregDerivs(int initStateId) {
    
    // get qureg ids (one for each var)
    int* quregIds;
    int numQuregs;
    WSGetInteger32List(stdlink, &quregIds, &numQuregs); // must free
    
    // get circuit indices of variables (ordered by return order)
    int* varOpInds;
    int numVars;
    WSGetInteger32List(stdlink, &varOpInds, &numVars); // must free

    // get ansatz circuit from MMA link
    int numOps;
    int *opcodes, *ctrls, *numCtrlsPerOp, *targs, *numTargsPerOp, *numParamsPerOp;
    qreal* params;
    int totalNumCtrls, totalNumTargs, totalNumParams; // these fields are only needed by clean-up
    local_loadCircuitFromMMA(
        &numOps, &opcodes, &ctrls, &numCtrlsPerOp, 
        &targs, &numTargsPerOp, &params, &numParamsPerOp,
        &totalNumCtrls, &totalNumTargs, &totalNumParams); // must free
    
    // get derivatives of any unitary matrices present in ansatz
    qreal* unitaryDerivs;
    int numElems;
    WSGetReal64List(stdlink, &unitaryDerivs, &numElems); // must free
    
    try {
        // validate inputs (note varOpInds is already validated by MMA caller)
        if (numQuregs != numVars)
            throw QuESTException("", "An equal number of quregs as variables must be passed.");
        
        // check initQureg is a created state-vector
        local_throwExcepIfQuregNotCreated(initStateId);
        int numQb = quregs[initStateId].numQubitsRepresented;
        
        /*
        if (quregs[initStateId].isDensityMatrix)
            throw QuESTException("", "Density matrices are not yet supported.");
        */
        
        for (int i=0; i < numQuregs; i++) {
            // check all quregs are created
            local_throwExcepIfQuregNotCreated(quregIds[i]);
            
            // check initStateId is not in derivatives
            if (quregIds[i] == initStateId)
                throw QuESTException("", "derivQuregs must not contain initQureg.");

            // check all quregs are the same size
            if (quregs[quregIds[i]].numQubitsRepresented != numQb)
                throw QuESTException("", "All qureg dimensions (of initQureg, and derivQuregs) must match.");
                
            // check all quregs are state-vectors
            /*
            if (quregs[quregIds[i]].isDensityMatrix)
                throw QuESTException("", "Density matrices are not yet supported.");
            */
        }
            
        // set all quregs to the initial state
        for (int q=0; q < numQuregs; q++)
            cloneQureg(quregs[quregIds[q]], quregs[initStateId]); // throw precluded by above validation
            
        local_getDerivativeQuregs(
            quregIds, varOpInds, numVars, numOps, opcodes, ctrls, numCtrlsPerOp, 
            targs, numTargsPerOp, params, numParamsPerOp, unitaryDerivs); // throws
            
        // return
        WSPutInteger(stdlink, initStateId);
        
        // proceed to clean-up afer catch
        
    } catch (QuESTException& err) {
        
        // report error, depending on type
        if (err.thrower == "")
            local_sendErrorAndFail("CalcQuregDerivs", err.message);
        else if (err.thrower == "Abort")
            local_sendErrorAndAbort("CalcQuregDerivs", err.message);
        else 
            local_sendErrorAndFail("CalcQuregDerivs", 
                "Problem in " + err.thrower + ": " + err.message);
                
        // proceed to clean-up below
    }
    
    // clean-up, even if errors have been sent to MMA
    WSReleaseInteger32List(stdlink, quregIds, numQuregs);
    WSReleaseInteger32List(stdlink, varOpInds, numVars);
    WSReleaseReal64List(stdlink, unitaryDerivs, numElems);
    local_freeCircuit(
        opcodes, ctrls, numCtrlsPerOp, targs, 
        numTargsPerOp, params, numParamsPerOp,
        numOps, totalNumCtrls, totalNumTargs, totalNumParams);
}

