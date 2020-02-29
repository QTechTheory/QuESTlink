/* @file quest_link.c
 * This file contains the C code for bridging QuEST to Mathematica.
 * While this file contains most of the necessary C code, there are some 
 * minor changes to the QuEST backend - overriding QuEST_validation's 
 * exitWithError, removing some validation, and adding an isCreated field 
 * to the Qureg struct.
 * 
 * The functions defined in this file are grouped as either 'local', 
 * 'callable', 'wrapper', 'internal':
 * - 'local' functions are called from only inside this file.
 * - 'callable' functions are directly callable by a user in MMA.
 * - 'wrapper' functions wrap a specific QuEST function and are directly 
 *    callable by a user in MMA.
 * - 'internal' functions are callable by the MMA, but shouldn't be called
 *    directly by a user - instead, they are wrapped by MMA-defined functions
 * 
 * Note the actual function names called from the MMA code differ from those 
 * below, as specified in quest_templates.tm
 *
 * There are a few layers of error-codes used in this file.
 * callable_ and wrapper_ functions often return -1 directly to the user,
 * to indicate user-error.
 * local_ functions return 0 to indicate an error, the explanation of which has 
 * been written to errorMsgBuffer.
 * internal_ functions with a Manual return must raise the error and push 
 * some kind of failed output (Abort, $Failed, -1) and avoid also returning a 
 * C primitive which will choke the pipeline.
 */


#include "wstp.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <QuEST.h>
#include <string>
#include <vector>

/*
 * PI constant needed for (multiControlled) sGate and tGate
 */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Codes for specifiying circuits from MMA
 */
#define OPCODE_H 0
#define OPCODE_X 1
#define OPCODE_Y 2
#define OPCODE_Z 3
#define OPCODE_Rx 4
#define OPCODE_Ry 5
#define OPCODE_Rz 6
#define OPCODE_R 7
#define OPCODE_S 8
#define OPCODE_T 9
#define OPCODE_U 10
#define OPCODE_Deph 11
#define OPCODE_Depol 12
#define OPCODE_Damp 13
#define OPCODE_SWAP 14
#define OPCODE_M 15
#define OPCODE_P 16
#define OPCODE_Kraus 17
#define OPCODE_G 18

/**
 * The initial number of maximum quregs represnetable before qureg-list resizing
 */
#define INIT_MAX_NUM_QUREGS 1000

/**
 * Max number of target and control qubits which can be specified 
 * for an individual gate 
 */
#define MAX_NUM_TARGS_CTRLS 100

/**
 * Global instance of QuESTEnv, created when MMA is linked.
 */
QuESTEnv env;

/**
 * Collection of instantiated Quregs
 */
std::vector<Qureg> quregs;
std::vector<bool> quregIsCreated;

/**
 * Buffer for creating error messages
 */
static char errorMsgBuffer[1000];

/**
 * Reports an error message to MMA without aborting
 */
void local_sendErrorToMMA(std::string err_msg) {
    WSPutFunction(stdlink, "EvaluatePacket", 1);
    WSPutFunction(stdlink, "Echo", 2);
    WSPutString(stdlink, err_msg.c_str());
    WSPutString(stdlink, "Error: ");
    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
}

void local_sendErrorMsgBufferToMMA(void) {
    local_sendErrorToMMA(errorMsgBuffer);
}

void local_sendQuregNotCreatedError(int id) {
    sprintf(errorMsgBuffer, "qureg (with id %d) has not been created.", id);
    local_sendErrorMsgBufferToMMA();
}

int local_writeToErrorMsgBuffer(std::string msg, ...) {
    va_list argp;
    va_start(argp, msg);
    vsprintf(errorMsgBuffer, msg.c_str(), argp);
    va_end(argp);
    return 0;
}

bool local_isQuregCreated(int id) {
    if (id < 0 || id >= (int) quregs.size())
        return false;
    return local_isQuregCreated(id);
}

size_t local_getNextQuregID(void) {
    size_t id;
    
    // check for next id
    for (id=0; id < quregs.size(); id++)
        if (!quregIsCreated[id])
            return id;
            
    // if none are available, make more space (using an empty Qureg)
    Qureg blank;
    id = quregs.size();
    quregs.push_back(blank);
    quregIsCreated.push_back(false);
    return id;
}

void wrapper_createQureg(int numQubits) {
    size_t id = local_getNextQuregID();
    quregs[id] = createQureg(numQubits, env);
    quregIsCreated[id] = true;
    WSPutInteger(stdlink, id);
}
void wrapper_createDensityQureg(int numQubits) {
    size_t id = local_getNextQuregID();
    quregs[id] = createDensityQureg(numQubits, env);
    quregIsCreated[id] = true;
    WSPutInteger(stdlink, id);
}
void wrapper_destroyQureg(int id) {
    if (local_isQuregCreated(id)) {
        destroyQureg(quregs[id], env);
        quregIsCreated[id] = false;
        WSPutInteger(stdlink, id);
    } else
        local_sendQuregNotCreatedError(id);
}

void callable_createQuregs(int numQubits, int numQuregs) {
    int ids[numQuregs];
    for (int i=0; i < numQuregs; i++) {
        int id = local_getNextQuregID();
        ids[i] = id;
        quregs[id] = createQureg(numQubits, env);
        quregIsCreated[id] = true;
    }
    WSPutIntegerList(stdlink, ids, numQuregs);
}

void callable_createDensityQuregs(int numQubits, int numQuregs) {
    int ids[numQuregs];
    for (int i=0; i < numQuregs; i++) {
        int id = local_getNextQuregID();
        ids[i] = id;
        quregs[id] = createDensityQureg(numQubits, env);
        quregIsCreated[id] = true;
    }
    WSPutIntegerList(stdlink, ids, numQuregs);
} 



/** getters */
void internal_getAmp(int quregID) {
    Qureg qureg = quregs[quregID];
    if (!local_isQuregCreated(quregID)) {
        local_sendQuregNotCreatedError(quregID);
        WSPutSymbol(stdlink, "$Failed");
        return;
    }
        
    // get args
    long int row;
    long int col;
    WSGetInteger64(stdlink, &row);
    WSGetInteger64(stdlink, &col);
    
    // ensure user supplied the correct number of args for the qureg type
    if (qureg.isDensityMatrix && col==-1) {
        local_sendErrorToMMA("GetAmp was called on a density matrix without suppling row and column.");
        WSPutSymbol(stdlink, "$Failed");
        return;
    }
    if (!qureg.isDensityMatrix && col!=-1) {
        local_sendErrorToMMA("GetAmp was called on a state-vector, yet a column index was supplied.");
        WSPutSymbol(stdlink, "$Failed");
        return;
    }
    
    // fetch amp (can throw internal QuEST errors)
    Complex amp;
    if (qureg.isDensityMatrix)
        amp = getDensityAmp(qureg, row, col);
    else
        amp = getAmp(qureg, row);
    
    // return as complex number
    WSPutFunction(stdlink, "Complex", 2);
    WSPutReal64(stdlink, amp.real);
    WSPutReal64(stdlink, amp.imag);
}

void callable_isDensityMatrix(int quregID) {
    Qureg qureg = quregs[quregID];
    if (!local_isQuregCreated(quregID)) {
        local_sendQuregNotCreatedError(quregID);
        WSPutSymbol(stdlink, "$Failed");
        return;
    }
    WSPutInteger(stdlink, qureg.isDensityMatrix);
}




/** initial states */

void wrapper_initZeroState(int id) {
    if (local_isQuregCreated(id))
        initZeroState(quregs[id]);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_initPlusState(int id) {
    if (local_isQuregCreated(id))
        initPlusState(quregs[id]);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_initClassicalState(int id, int stateInd) {
    if (local_isQuregCreated(id))
        initClassicalState(quregs[id], stateInd);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_initPureState(int quregID, int pureID) {
    if (local_isQuregCreated(quregID) && local_isQuregCreated(pureID))
        initPureState(quregs[quregID], quregs[pureID]);
    else if (!local_isQuregCreated(quregID))
        local_sendQuregNotCreatedError(quregID);
    else// if (!local_isQuregCreated(pureID))
        local_sendQuregNotCreatedError(pureID);
    WSPutInteger(stdlink, quregID);
}
void wrapper_initStateFromAmps(int quregID, qreal* reals, long l1, qreal* imags, long l2) {
    Qureg qureg = quregs[quregID];
    
    if (!local_isQuregCreated(quregID))
        local_sendQuregNotCreatedError(quregID);
    else if (l1 != l2 || l1 != qureg.numAmpsTotal)
        local_sendErrorToMMA("incorrect number of amplitudes supplied! State has not been changed.");
    else
        initStateFromAmps(qureg, reals, imags);
    WSPutInteger(stdlink, quregID);
}
void wrapper_cloneQureg(int outID, int inID) {
    if (!local_isQuregCreated(outID))
        local_sendQuregNotCreatedError(outID);
    else if (!local_isQuregCreated(inID))
        local_sendQuregNotCreatedError(inID);
    else
        cloneQureg(quregs[outID], quregs[inID]);
    WSPutInteger(stdlink, outID);
}





/** noise */

void wrapper_mixDephasing(int id, int qb1, qreal prob) {
    if (local_isQuregCreated(id))
        mixDephasing(quregs[id], qb1, prob);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_mixTwoQubitDephasing(int id, int qb1, int qb2, qreal prob) {
    if (local_isQuregCreated(id))
        mixTwoQubitDephasing(quregs[id], qb1, qb2, prob);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_mixDepolarising(int id, int qb1, qreal prob) {
    if (local_isQuregCreated(id))
        mixDepolarising(quregs[id], qb1, prob);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_mixTwoQubitDepolarising(int id, int qb1, int qb2, qreal prob) {
    if (local_isQuregCreated(id))
        mixTwoQubitDepolarising(quregs[id], qb1, qb2, prob);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}
void wrapper_mixDamping(int id, int qb, qreal prob) {
    if (local_isQuregCreated(id))
        mixDamping(quregs[id], qb, prob);
    else
        local_sendQuregNotCreatedError(id);
    WSPutInteger(stdlink, id);
}


/* calculations */

void wrapper_calcProbOfOutcome(int id, int qb, int outcome) {
    if (local_isQuregCreated(id)) {
        WSPutReal64(stdlink, calcProbOfOutcome(quregs[id], qb, outcome));
        return;
    }
    else {
        local_sendQuregNotCreatedError(id);
        WSPutInteger(stdlink, -1);
        return;
    }
}
void wrapper_calcFidelity(int id1, int id2) {
    if (!quregIsCreated[id1]) {
        local_sendQuregNotCreatedError(id1);
        WSPutInteger(stdlink, -1);
        return;
    }
    if (!quregIsCreated[id2]) {
        local_sendQuregNotCreatedError(id2);
        WSPutInteger(stdlink, -1);
        return;
    }
    WSPutReal64(stdlink, calcFidelity(quregs[id1], quregs[id2]));
}
void wrapper_calcInnerProduct(int id1, int id2) {
    if (!quregIsCreated[id1]) {
        local_sendQuregNotCreatedError(id1);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    if (!quregIsCreated[id2]) {
        local_sendQuregNotCreatedError(id2);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
        
    Complex res = calcInnerProduct(quregs[id1], quregs[id2]);
    WSPutFunction(stdlink, "Complex", 2);
    WSPutReal64(stdlink, res.real);
    WSPutReal64(stdlink, res.imag);
}
void wrapper_calcDensityInnerProduct(int id1, int id2) {
    if (!quregIsCreated[id1]) {
        local_sendQuregNotCreatedError(id1);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    if (!quregIsCreated[id2]) {
        local_sendQuregNotCreatedError(id2);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    qreal res = calcDensityInnerProduct(quregs[id1], quregs[id2]);
    WSPutReal64(stdlink, res);
}
void wrapper_calcPurity(int id) {
    if (!local_isQuregCreated(id)) {
        local_sendQuregNotCreatedError(id);
        WSPutInteger(stdlink, -1);
        return;
    }
    WSPutReal64(stdlink, calcPurity(quregs[id]));
}
void wrapper_calcTotalProb(int id) {
    if (!local_isQuregCreated(id)) {
        local_sendQuregNotCreatedError(id);
        WSPutInteger(stdlink, -1);
        return;
    }
    WSPutReal64(stdlink, calcTotalProb(quregs[id]));
}
void wrapper_calcHilbertSchmidtDistance(int id1, int id2) {
    if (!quregIsCreated[id1]) {
        local_sendQuregNotCreatedError(id1);
        WSPutInteger(stdlink, -1);
        return;
    }
    if (!quregIsCreated[id2]) {
        local_sendQuregNotCreatedError(id2);
        WSPutInteger(stdlink, -1);
        return;
    }
    WSPutReal64(stdlink, calcHilbertSchmidtDistance(quregs[id1], quregs[id2]));
}


/* other modifications */

void wrapper_collapseToOutcome(int id, int qb, int outcome) {
    if (local_isQuregCreated(id)) {
        collapseToOutcome(quregs[id], qb, outcome);
        WSPutInteger(stdlink, id);
    } else {
        local_sendQuregNotCreatedError(id);
        WSPutInteger(stdlink, -1);
    }
}


/* circuit execution */

int local_gateUnsupportedError(std::string gate) {
    return local_writeToErrorMsgBuffer("the gate '%s' is not supported.", gate);
}

int local_gateWrongNumParamsError(std::string gate, int wrongNumParams, int rightNumParams) {
    return local_writeToErrorMsgBuffer(
        "the gate '%s' accepts %d parameters, but %d were passed.",
        gate, rightNumParams, wrongNumParams);
}

/* rightNumTargs is a string so that it can be multiple e.g. "1 or 2" */
int local_gateWrongNumTargsError(std::string gate, int wrongNumTargs, std::string rightNumTargs) {
    return local_writeToErrorMsgBuffer(
        "the gate '%s' accepts %s, but %d were passed.",
        gate, rightNumTargs, wrongNumTargs);
}

qreal local_getMaxValidNoiseProb(int opcode, int numQubits) {
    if (opcode == OPCODE_Damp)
        return 1.0;
    if (opcode == OPCODE_Deph) {
        if (numQubits == 1)
            return 1.0/2.0;
        if (numQubits == 2)
            return 3.0/4.0;
    }
    if (opcode == OPCODE_Depol) {
        if (numQubits == 1)
            return 3.0/4.0;
        if (numQubits == 2)
            return 15.0/16.0;
    }
    return -1; // should never be reached
}

int local_isValidProb(int opcode, int numQubits, qreal prob) {
    return (prob > 0 && prob < local_getMaxValidNoiseProb(opcode, numQubits));
}

int local_noiseInvalidProbError(int opcode, int numQubits, qreal prob) {
    std::string opStr = "";
    if (opcode == OPCODE_Deph) opStr = "dephasing";
    if (opcode == OPCODE_Depol) opStr = "depolarising";
    if (opcode == OPCODE_Damp) opStr = "amplitude damping";
    qreal maxProb = local_getMaxValidNoiseProb(opcode, numQubits);
    
    return local_writeToErrorMsgBuffer(
        "%d-qubit %s was applied with probability %g which is outside its accepted range of [0, %g].", 
        numQubits, opStr, prob, maxProb);
}

int* local_prepareCtrlCache(int* ctrls, int ctrlInd, int numCtrls, int addTarg) {
    static int ctrlCache[MAX_NUM_TARGS_CTRLS]; 
    for (int i=0; i < numCtrls; i++)
        ctrlCache[i] = ctrls[ctrlInd + i];
    if (addTarg != -1)
        ctrlCache[numCtrls] = addTarg;
    return ctrlCache;
}

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



/* returns 1 if successful, else 0, upon which caller must error and abort.
 * if returns 0, errors messages will be written to global errorMsgBuffer
 * @param mesOutcomeCache may be NULL
 * @param finalCtrlInd, finalTargInd and finalParamInd are modified to point to 
 *  the final values of ctrlInd, targInd and paramInd, after the #numOps operation 
 *  has been applied. If #numOps isn't smaller than the actual length of the 
 *  circuit arrays, these indices will point out of bounds.
 */
int local_applyGates(
    Qureg qureg, 
    int numOps, int* opcodes, 
    int* ctrls, int* numCtrlsPerOp, 
    int* targs, int* numTargsPerOp, 
    qreal* params, int* numParamsPerOp,
    int* mesOutcomeCache,
    int* finalCtrlInd, int* finalTargInd, int* finalParamInd
    ) {
        
    int ctrlInd = 0;
    int targInd = 0;
    int paramInd = 0;
    int mesInd = 0;
    
    // attempt to apply each gate
    for (int opInd=0; opInd < numOps; opInd++) {
        
        //@TODO: benchmark time cost of abort checking!
        
        // check whether the user has tried to abort
        if (WSMessageReady(stdlink)) {
            int code, arg;
            WSGetMessage(stdlink, &code, &arg);
            if (code == WSTerminateMessage || code == WSInterruptMessage || 
                code == WSAbortMessage     || code == WSImDyingMessage) {
                return local_writeToErrorMsgBuffer("Circuit simulation aborted.");
            }
        }

        // get gate info
        int op = opcodes[opInd];
        int numCtrls = numCtrlsPerOp[opInd];
        int numTargs = numTargsPerOp[opInd];
        int numParams = numParamsPerOp[opInd];
        
        switch(op) {
            
            case OPCODE_H :
                if (numParams != 0)
                    return local_gateWrongNumParamsError("Hadamard", numParams, 0);
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled Hadamard");
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("Hadamard", numTargs, "1 target");
                hadamard(qureg, targs[targInd]);
                break;
                
            case OPCODE_S :
                if (numParams != 0)
                    return local_gateWrongNumParamsError("S gate", numParams, 0);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("S gate", numTargs, "1 target");
                if (numCtrls == 0)
                    sGate(qureg, targs[targInd]);
                else {
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledPhaseShift(qureg, ctrlCache, numCtrls+1, M_PI/2);
                }
                break;
                
            case OPCODE_T : {
                if (numParams != 0)
                    return local_gateWrongNumParamsError("T gate", numParams, 0);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("T gate", numTargs, "1 target");
                if (numCtrls == 0)
                    tGate(qureg, targs[targInd]);
                else {
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledPhaseShift(qureg, ctrlCache, numCtrls+1, M_PI/4);
                }
            }
                break;
        
            case OPCODE_X : {
                if (numParams != 0)
                    return local_gateWrongNumParamsError("X", numParams, 0);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("X", numTargs, "1 target");
                if (numCtrls == 0)
                    pauliX(qureg, targs[targInd]);
                else if (numCtrls == 1)
                    controlledNot(qureg, ctrls[ctrlInd], targs[targInd]);
                else {
                    ComplexMatrix2 u = {
                        .real={{0,1},{1,0}},
                        .imag={{0}}};
                    multiControlledUnitary(qureg, &ctrls[ctrlInd], numCtrls, targs[targInd], u);
                }
            }
                break;
                
            case OPCODE_Y :
                if (numParams != 0)
                    return local_gateWrongNumParamsError("Y", numParams, 0);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("Y", numTargs, "1 target");
                if (numCtrls == 0)
                    pauliY(qureg, targs[targInd]);
                else if (numCtrls == 1)
                    controlledPauliY(qureg, ctrls[ctrlInd], targs[targInd]);
                else
                    return local_gateUnsupportedError("controlled Y");
                break;
                
            case OPCODE_Z : {
                if (numParams != 0)
                    return local_gateWrongNumParamsError("Z", numParams, 0);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("Z", numTargs, "1 target");
                if (numCtrls == 0)
                    pauliZ(qureg, targs[targInd]);
                else {
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledPhaseFlip(qureg, ctrlCache, numCtrls+1);
                }
            }
                break;
        
            case OPCODE_Rx :
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Rx", numParams, 1);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("Rx", numTargs, "1 target");
                if (numCtrls == 0)
                    rotateX(qureg, targs[targInd], params[paramInd]);
                else if (numCtrls == 1)
                    controlledRotateX(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]);
                else
                    return local_gateUnsupportedError("multi-controlled Rotate X");
                break;
                
            case OPCODE_Ry :
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Ry", numParams, 1);
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("Ry", numTargs, "1 target");
                if (numCtrls == 0)
                    rotateY(qureg, targs[targInd], params[paramInd]);
                else if (numCtrls == 1)
                    controlledRotateY(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]);
                else
                    return local_gateUnsupportedError("multi-controlled Rotate Y");
                break;
                
            case OPCODE_Rz :
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Rz", numParams, 1);
                if (numCtrls > 1)
                    return local_gateUnsupportedError("multi-controlled Rotate Z");
                if (numCtrls == 1 && numTargs > 1)
                    return local_gateUnsupportedError("multi-controlled multi-rotateZ");
                if (numTargs == 1) {
                    if (numCtrls == 0)
                        rotateZ(qureg, targs[targInd], params[paramInd]);
                    if (numCtrls == 1)
                        controlledRotateZ(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]);
                } else
                    multiRotateZ(qureg, &targs[targInd], numTargs, params[paramInd]);
                break;
                
            case OPCODE_R: {
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled multi-rotate-Pauli");
                if (numTargs != numParams-1) {
                    return local_writeToErrorMsgBuffer(
                        "An internel error in R occured! "
                        "The quest_link received an unequal number of Pauli codes (%d) and target qubits! (%d)",
                        numParams-1, numTargs);
                }
                enum pauliOpType paulis[MAX_NUM_TARGS_CTRLS]; 
                for (int p=0; p < numTargs; p++)
                    paulis[p] = (pauliOpType) ((int) params[paramInd+1+p]);
                multiRotatePauli(qureg, &targs[targInd], paulis, numTargs, params[paramInd]);
            }
                break;
            
            case OPCODE_U : {
                if (numTargs == 1 && numParams != 2*2*2)
                    return local_writeToErrorMsgBuffer("single qubit U accepts only 2x2 matrices");
                if (numTargs == 2 && numParams != 4*4*2)
                    return local_writeToErrorMsgBuffer("two qubit U accepts only 4x4 matrices");
                if (numTargs != 1 && numTargs != 2)
                    return local_gateWrongNumTargsError("U", numTargs, "1 or 2 targets");
                
                if (numTargs == 1) {
                    ComplexMatrix2 u = local_getMatrix2FromFlatList(&params[paramInd]);
                    if (numCtrls == 0)
                        unitary(qureg, targs[targInd], u);
                    else
                        multiControlledUnitary(qureg, &ctrls[ctrlInd], numCtrls, targs[targInd], u);
                }
                else if (numTargs == 2) {
                    ComplexMatrix4 u = local_getMatrix4FromFlatList(&params[paramInd]);
                    if (numCtrls == 0)
                        twoQubitUnitary(qureg, targs[targInd], targs[targInd+1], u);
                    else
                        multiControlledTwoQubitUnitary(qureg, &ctrls[ctrlInd], numCtrls, targs[targInd], targs[targInd+1], u);
                }
            }
                break;
                
            case OPCODE_Deph :
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Dephasing", numParams, 1);
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled dephasing");
                if (numTargs != 1 && numTargs != 2)
                    return local_gateWrongNumTargsError("Dephasing", numTargs, "1 or 2 targets");
                if (params[paramInd] == 0)
                    break;
                if (!local_isValidProb(op, numTargs, params[paramInd]))
                    return local_noiseInvalidProbError(op, numTargs, params[paramInd]);
                if (numTargs == 1)
                    mixDephasing(qureg, targs[targInd], params[paramInd]);
                if (numTargs == 2)
                    mixTwoQubitDephasing(qureg, targs[targInd], targs[targInd+1], params[paramInd]);
                break;
                
            case OPCODE_Depol :
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Depolarising", numParams, 1);
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled depolarising");
                if (numTargs != 1 && numTargs != 2)
                    return local_gateWrongNumTargsError("Depolarising", numTargs, "1 or 2 targets");
                if (params[paramInd] == 0)
                    break;
                if (!local_isValidProb(op, numTargs, params[paramInd]))
                    return local_noiseInvalidProbError(op, numTargs, params[paramInd]);
                if (numTargs == 1)
                    mixDepolarising(qureg, targs[targInd], params[paramInd]);
                if (numTargs == 2)
                    mixTwoQubitDepolarising(qureg, targs[targInd], targs[targInd+1], params[paramInd]);
                break;
                
            case OPCODE_Damp :
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Damping", numParams, 1);
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled damping");
                if (numTargs != 1)
                    return local_gateWrongNumTargsError("Damping", numTargs, "1 target");
                if (params[paramInd] == 0)
                    break;
                if (!local_isValidProb(op, numTargs, params[paramInd]))
                    return local_noiseInvalidProbError(op, numTargs, params[paramInd]);
                mixDamping(qureg, targs[targInd], params[paramInd]);
                break;
                
            case OPCODE_SWAP: {
                if (numParams != 0)
                    return local_gateWrongNumParamsError("SWAP", numParams, 0);
                if (numTargs != 2)
                    return local_gateWrongNumTargsError("Depolarising", numTargs, "2 targets");
                if (numCtrls == 0) {
                    controlledNot(qureg, targs[targInd],   targs[targInd+1]);
                    controlledNot(qureg, targs[targInd+1], targs[targInd  ]);
                    controlledNot(qureg, targs[targInd],   targs[targInd+1]);
                } else {
                    ComplexMatrix2 u = {
                        .real={{0,1},{1,0}},
                        .imag={{0}}};
                    int* ctrlCache = local_prepareCtrlCache(ctrls, ctrlInd, numCtrls, targs[targInd]);
                    multiControlledUnitary(qureg, ctrlCache, numCtrls+1, targs[targInd+1], u);
                    ctrlCache[numCtrls] = targs[targInd+1];
                    multiControlledUnitary(qureg, ctrlCache, numCtrls+1, targs[targInd], u);
                    ctrlCache[numCtrls] = targs[targInd];
                    multiControlledUnitary(qureg, ctrlCache, numCtrls+1, targs[targInd+1], u);
                }
            }
                break;
                
            case OPCODE_M: {
                if (numParams != 0)
                    return local_gateWrongNumParamsError("M", numParams, 0);
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled measurement");
                for (int q=0; q < numTargs; q++) {
                    int outcomeVal = measure(qureg, targs[targInd+q]);
                    if (mesOutcomeCache != NULL)
                        mesOutcomeCache[mesInd++] = outcomeVal;
                }
            }
                break;
            
            case OPCODE_P:
                if (numParams != 1 && numParams != numTargs) {
                    return local_writeToErrorMsgBuffer(
                        "P[outcomes] specified a different number of binary outcomes (%d) than target qubits (%d)!",
                        numParams, numTargs);
                }
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled projector");
                if (numParams > 1)
                    for (int q=0; q < numParams; q++)
                        collapseToOutcome(qureg, targs[targInd+q], (int) params[paramInd+q]);
                else {
                    // check value isn't impossibly high
                    if (params[paramInd] >= (1LL << numTargs))
                        return local_writeToErrorMsgBuffer(
                            "P[%d] was applied to %d qubits and exceeds their maximum represented value of %lld.",
                            (int) params[paramInd], numTargs, (1LL << numTargs));
                    // work out each bit outcome and apply; right most (least significant) bit acts on right-most target
                    for (int q=0; q < numTargs; q++)
                        collapseToOutcome(qureg, targs[targInd+numTargs-q-1], (((int) params[paramInd]) >> q) & 1);
                }
                break;
                
            case OPCODE_Kraus: {
                ; // empty post-label statement, courtesy of weird C99 standard
                int numKrausOps = (int) params[paramInd];
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled Kraus map");
                if (numTargs != 1 && numTargs != 2)
                    return local_gateWrongNumTargsError("Kraus map", numTargs, "1 or 2 targets");
                if ((numKrausOps < 1) ||
                    (numTargs == 1 && numKrausOps > 4 ) ||
                    (numTargs == 2 && numKrausOps > 16))
                    return local_writeToErrorMsgBuffer(
                        "%d operators were passed to single-qubit Krauss[ops], which accepts only >0 and <=%d operators!",
                        numKrausOps, (numTargs==1)? 4:16);
                if (numTargs == 1 && (numParams-1) != 2*2*2*numKrausOps)
                    return local_writeToErrorMsgBuffer("one-qubit Kraus expects 2-by-2 matrices!");
                if (numTargs == 2 && (numParams-1) != 4*4*2*numKrausOps)
                    return local_writeToErrorMsgBuffer("two-qubit Kraus expects 4-by-4 matrices!");

                if (numTargs == 1) {
                    ComplexMatrix2 krausOps[4];
                    int opElemInd = 1 + paramInd;
                    for (int n=0; n < numKrausOps; n++)
                        krausOps[n] = local_getMatrix2FromFlatList(&params[opElemInd + 2*2*2*n]);
                    mixKrausMap(qureg, targs[targInd], krausOps, numKrausOps);
                } 
                else if (numTargs == 2) {
                    ComplexMatrix4 krausOps[16];
                    int opElemInd = 1 + paramInd;
                    for (int n=0; n < numKrausOps; n++)
                        krausOps[n] = local_getMatrix4FromFlatList(&params[opElemInd + 2*4*4*n]);
                    mixTwoQubitKrausMap(qureg, targs[targInd], targs[targInd+1], krausOps, numKrausOps);
                }
            }
                break;
                
            case OPCODE_G : {
                if (numParams != 1)
                    return local_gateWrongNumParamsError("Global phase", numParams, 1);
                if (numCtrls != 0)
                    return local_gateUnsupportedError("controlled global phase");
                if (numTargs != 0)
                    return local_gateWrongNumTargsError("Global phase", numTargs, "0 targets");
                if (params[paramInd] == 0)
                    break;
                // phase does not change density matrices
                if (!qureg.isDensityMatrix) {
                    Complex zero = (Complex) {.real=0, .imag=0};
                    Complex fac = (Complex) {.real=cos(params[paramInd]), .imag=sin(params[paramInd])};
                    setWeightedQureg(zero, qureg, zero, qureg, fac, qureg); // exp(i param)|qureg>
                }
            }
                break;
                
            default:            
                return local_writeToErrorMsgBuffer("circuit contained an unknown gate.");
        }
        ctrlInd += numCtrls;
        targInd += numTargs;
        paramInd += numParams;
    }
    
    // update final pointers
    *finalCtrlInd = ctrlInd;
    *finalTargInd = targInd;
    *finalParamInd = paramInd;
    
    // indicate success 
    return 1;
}

void local_loadCircuitFromMMA(
    int* numOps, int** opcodes, int** ctrls, int** numCtrlsPerOp, 
    int** targs, int** numTargsPerOp, qreal** params, int** numParamsPerOp) {

    int totalNumCtrls, totalNumTargs, totalNumParams;

    WSGetInteger32List(stdlink, opcodes, numOps);
    
    WSGetInteger32List(stdlink, ctrls, &totalNumCtrls);
    WSGetInteger32List(stdlink, numCtrlsPerOp, numOps);

    WSGetInteger32List(stdlink, targs, &totalNumTargs);
    WSGetInteger32List(stdlink, numTargsPerOp, numOps);

    WSGetReal64List(stdlink, params, &totalNumParams);
    WSGetInteger32List(stdlink, numParamsPerOp, numOps);
}

/** 
 * Applies a given circuit to the identified qureg.
 * The circuit is expressed as lists of opcodes (identifying gates),
 * the total flat sequence control qubits, a list denoting how many of
 * the control qubits apply to each operation, their target qubits (flat list),
 * a list denotating how many targets each operation has, their parameters 
 * (flat list) and a list denoting how many params each operation has.
 * Returns a list of measurement outcome of any performed measurements in the circuit.
 * The original qureg of the state is restored when this function
 * is aborted by the calling MMA, or aborted due to encountering
 * an invalid gate. In this case, Abort[] is returned.
 * However, a user error caught by the QuEST backend
 * (e.g. same target and control qubit) will result in the link being
 * destroyed.
 */
void internal_applyCircuit(int id) {
    
    // get arguments from MMA link
    int numOps;
    int *opcodes, *ctrls, *numCtrlsPerOp, *targs, *numTargsPerOp, *numParamsPerOp;
    qreal* params;
    local_loadCircuitFromMMA(&numOps, &opcodes, &ctrls, &numCtrlsPerOp, &targs, &numTargsPerOp, &params, &numParamsPerOp);
        
    // ensure qureg exists
    if (!local_isQuregCreated(id)) {
        local_sendQuregNotCreatedError(id);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    
    Qureg qureg = quregs[id];
    
    // backup of initial state in case of abort
    Qureg backup = createCloneQureg(qureg, env);
    
    // count the total number of measurements performed in a circuit
    int totalNumMesGates = 0;
    int totalNumMeasurements = 0;
    for (int opInd=0; opInd < numOps; opInd++)
        if (opcodes[opInd] == OPCODE_M) {
            totalNumMesGates++;
            totalNumMeasurements += numTargsPerOp[opInd];
        }
        
    // prepare records of measurement outcomes
    int* mesOutcomeCache = (int*) malloc(totalNumMeasurements * sizeof(int));
    int mesInd = 0;
    
    // apply the circuit
    int finalCtrlInd, finalTargInd, finalParamInd; // ignored final inds
    int success = local_applyGates(
        qureg, numOps, opcodes, 
        ctrls, numCtrlsPerOp, 
        targs, numTargsPerOp, 
        params, numParamsPerOp,
        mesOutcomeCache,
        &finalCtrlInd, &finalTargInd, &finalParamInd);
    
    // if circuit contained no errors...
    if (success) {
        // return lists of measurement outcomes
        mesInd = 0;
        WSPutFunction(stdlink, "List", totalNumMesGates);
        for (int opInd=0; opInd < numOps; opInd++) {
            if (opcodes[opInd] == OPCODE_M) {
                WSPutFunction(stdlink, "List", numTargsPerOp[opInd]);
                for (int i=0; i < numTargsPerOp[opInd]; i++)
                    WSPutInteger(stdlink, mesOutcomeCache[mesInd++]);
            }
        }
    } else {
        // otherwise restore the qureg's original state and issue errors and Abort in Mathematica
        cloneQureg(qureg, backup);
        sprintf(
            errorMsgBuffer+strlen(errorMsgBuffer), 
            " Aborting circuit and restore qureg (id %d) to its original state", id);
        local_sendErrorMsgBufferToMMA();
        WSPutFunction(stdlink, "Abort", 0);
    }

    // clear data structures
    free(mesOutcomeCache);
    destroyQureg(backup, env);
}

/* quregs must be prior initialised and cloned to the initial state of the circuit.
 * returns 0 if each circuit and derivative was performed successful, else returns 
 * 1 and writes the error message to errorMsgBuffer
 */
int local_getDerivativeQuregs(
    // variable (to be differentiated) info
    int* quregIds, int* varOpInds, int numVars, 
    // circuit info
    int numOps, int* opcodes, 
    int* ctrls, int* numCtrlsPerOp, 
    int* targs, int* numTargsPerOp, 
    qreal* params, int* numParamsPerOp,
    // derivative matrices of general unitary gates in circuit
    qreal* unitaryDerivs) {
        
    // index of the first (real) element of the next unitary derivative
    int unitaryDerivInd = 0;

    // don't record measurement outcomes
    int* mesOutcomes = NULL;
    
    for (int v=0; v<numVars; v++) {
        Qureg qureg = quregs[quregIds[v]];
        int varOp = varOpInds[v];
        
        // indices AFTER last gate applied by circuit
        int finalCtrlInd, finalTargInd, finalParamInd;
        
        // apply only gates up to and including the to-be-differentiated gate,
        // unless that gate is the general unitary
        int op = opcodes[varOp];
        int diffGateWasApplied = (op != OPCODE_U);
        int success = local_applyGates(
            qureg, (diffGateWasApplied)? varOp+1 : varOp, opcodes, 
            ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp,
            mesOutcomes, &finalCtrlInd, &finalTargInd, &finalParamInd);
        if (!success)
            return success; // errorMsgBuffer updated by local_applyGates

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
        
        // choices of re-normalisation
        Complex negHalfI = (Complex) {.real=0, .imag=-0.5};
        Complex posI = (Complex) {.real=0, .imag=1};
        Complex zero = (Complex) {.real=0, .imag=0};
        Complex one = (Complex) {.real=1, .imag=0};
        
        Complex normFac;
        
        // ignore control qubits and apply gate Paulis incurred by differentiation 
        switch(op) {
            case OPCODE_Rx:
                for (int t=0; t < numTargs; t++) // multi-target X may be possible later 
                    pauliX(qureg, targs[t+finalTargInd]); 
                normFac = negHalfI;
                break;
            case OPCODE_Ry:
                for (int t=0; t < numTargs; t++) // multi-target Y may be possible later 
                    pauliY(qureg, targs[t+finalTargInd]); 
                normFac = negHalfI;
                break;
            case OPCODE_Rz:
                for (int t=0; t < numTargs; t++)
                    pauliZ(qureg, targs[t+finalTargInd]); 
                normFac = negHalfI;
                break;
            case OPCODE_R:
                for (int t=0; t < numTargs; t++) {
                    int pauliCode = (int) params[t+finalParamInd+1];
                    if (pauliCode == 1) pauliX(qureg, targs[t+finalTargInd]);
                    if (pauliCode == 2) pauliY(qureg, targs[t+finalTargInd]);
                    if (pauliCode == 3) pauliZ(qureg, targs[t+finalTargInd]);
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
                    applyOneQubitMatrix(qureg, targs[finalTargInd], u2);
                } else if (numTargs == 2) {
                    ComplexMatrix4 u4 = local_getMatrix4FromFlatList(&unitaryDerivs[unitaryDerivInd]);
                    unitaryDerivInd += 2*4*4;
                    applyTwoQubitMatrix(qureg, targs[finalTargInd], targs[finalTargInd+1], u4);
                }
                else {
                    // create a non-dynamic ComplexMatrixN instance 
                    // general matrix N unpacking here; can do static
                    return local_writeToErrorMsgBuffer("multi-qubit U deriv not yet implemented!");
                }
                normFac = one;
                break;
            default:            
                return local_writeToErrorMsgBuffer("Only Rx, Ry, Rz, R, U and their controlled gates may be differentiated!");
        }
        
        // differentiate control qubits by forcing them to 1, without renormalising
        for (int c=0; c<numCtrls; c++)
            projectToOne(qureg, ctrls[finalCtrlInd]);
        
        // adjust normalisation
        setWeightedQureg(zero, qureg, zero, qureg, normFac, qureg);
        
        // wind forward inds to point to the next gate 
        finalCtrlInd += numCtrls;
        finalTargInd += numTargs;
        finalParamInd += numParams;

        // apply the remainder of the circuit
        success = local_applyGates(
            qureg, numOps-(varOp+1), &opcodes[varOp+1], 
            &ctrls[  finalCtrlInd], &numCtrlsPerOp[varOp+1], 
            &targs[  finalTargInd], &numTargsPerOp[varOp+1], 
            &params[finalParamInd], &numParamsPerOp[varOp+1],
            mesOutcomes, &finalCtrlInd, &finalTargInd, &finalParamInd);
        if (!success)
            return success; // errorMsgBuffer updated by local_applyGates
    }
    
    // indicate success
    return 1;
}

void internal_calcQuregDerivs(int initStateId) {
    
    // get qureg ids (one for each var)
    int* quregIds;
    int numQuregs;
    WSGetInteger32List(stdlink, &quregIds, &numQuregs);
    
    // get circuit indices of variables (ordered by return order)
    int* varOpInds;
    int numVars;
    WSGetInteger32List(stdlink, &varOpInds, &numVars);

    // get ansatz circuit from MMA link
    int numOps;
    int *opcodes, *ctrls, *numCtrlsPerOp, *targs, *numTargsPerOp, *numParamsPerOp;
    qreal* params;
    local_loadCircuitFromMMA(&numOps, &opcodes, &ctrls, &numCtrlsPerOp, &targs, &numTargsPerOp, &params, &numParamsPerOp);
    
    qreal* unitaryDerivs;
    int numElems;
    WSGetReal64List(stdlink, &unitaryDerivs, &numElems);
    
    // check MMA-loaded args are valid
    if (numQuregs != numVars) {
        local_sendErrorToMMA("An equal number of quregs as variables must be passed.");
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    if (!local_isQuregCreated(initStateId)) {
        local_sendQuregNotCreatedError(initStateId);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    for (int i=0; i < numQuregs; i++)
        if (!quregIsCreated[quregIds[i]]) {
            local_sendQuregNotCreatedError(quregIds[i]);
            WSPutSymbol(stdlink, "$Failed");;
            return;
        }
    // varOpInds validated by MMA caller
            
    // set all quregs to the initial state
    for (int q=0; q < numQuregs; q++)
        cloneQureg(quregs[quregIds[q]], quregs[initStateId]);
            
    // compute derivatives 
    int success = local_getDerivativeQuregs(
        quregIds, varOpInds, numVars, 
        numOps, opcodes, 
        ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, 
        unitaryDerivs);
        
    if (!success) {
        local_sendErrorMsgBufferToMMA();
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    
    // need to send anything to fulfill MMA return
    WSPutInteger(stdlink, initStateId);
}

/* returns a real vector with ith element calcDensityInnerProduct(qureg[rhoId], qureg[omegaIds[i]]) */
void internal_calcDensityInnerProductsVector(int rhoId, int omegaIds[], long numOmegas) {
    
    // check quregs are instantiated
    if (!local_isQuregCreated(rhoId)) {
        local_sendQuregNotCreatedError(rhoId);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    for (int i=0; i<numOmegas; i++) {
        if (!quregIsCreated[omegaIds[i]]) {
            local_sendQuregNotCreatedError(omegaIds[i]);
            WSPutSymbol(stdlink, "$Failed");;
            return;
        }
    }
    
    // calculate inner products 
    qreal* prods = (qreal*) malloc(numOmegas * sizeof *prods);
    for (int i=0; i<numOmegas; i++)
        prods[i] = calcDensityInnerProduct(quregs[rhoId], quregs[omegaIds[i]]);

    // send result to MMA
    WSPutReal64List(stdlink, prods, numOmegas);
    
    free(prods);
}

/* returns vector with ith element <qureg[braId]|qureg[ketIds[i]]> */
void internal_calcInnerProductsVector(int braId, int ketIds[], long numKets) {
    
    // check quregs are instantiated
    if (!local_isQuregCreated(braId)) {
        local_sendQuregNotCreatedError(braId);
        WSPutSymbol(stdlink, "$Failed");;
        return;
    }
    for (int i=0; i<numKets; i++) {
        if (!quregIsCreated[ketIds[i]]) {
            local_sendQuregNotCreatedError(ketIds[i]);
            WSPutSymbol(stdlink, "$Failed");;
            return;
        }
    }
    
    // calculate inner products 
    qreal* vecRe = (qreal*) malloc(numKets * sizeof *vecRe);
    qreal* vecIm = (qreal*) malloc(numKets * sizeof *vecIm);
    for (int i=0; i<numKets; i++) {
        Complex val = calcInnerProduct(quregs[braId], quregs[ketIds[i]]);
        vecRe[i] = val.real;
        vecIm[i] = val.imag;
    }
    
    // send result to MMA
    WSPutFunction(stdlink, "List", 2);
    WSPutReal64List(stdlink, vecRe, numKets);
    WSPutReal64List(stdlink, vecIm, numKets);
    
    free(vecRe);
    free(vecIm);
}

/* returns real, symmetric matrix with ith jth element calcDensityInnerProduct(quregs[quregIds[i]], qureg[qurgIds[j]]) */
void internal_calcDensityInnerProductsMatrix(int quregIds[], long numQuregs) {
    
    // check all quregs are created
    for (int i=0; i<numQuregs; i++) {
        if (!quregIsCreated[quregIds[i]]) {
            local_sendQuregNotCreatedError(quregIds[i]);
            WSPutSymbol(stdlink, "$Failed");;
            return;
        }
    }
    
    // store real matrix as `nested pointers`
    long len = numQuregs * numQuregs;
    qreal* matr = (qreal*) malloc(len * sizeof *matr);
    
    for (int r=0; r<numQuregs; r++) {
        for (int c=0; c<numQuregs; c++) {
            if (c >= r) {
                qreal val = calcDensityInnerProduct(quregs[quregIds[r]], quregs[quregIds[c]]);
                matr[r*numQuregs + c] = val;
            } else {
                matr[r*numQuregs + c] = matr[c*numQuregs + r];
            }
        }
    }
    
    WSPutReal64List(stdlink, matr, len);
    free(matr);
}

/* returns Hermitian matrix with ith jth element <qureg[i]|qureg[j]> */
void internal_calcInnerProductsMatrix(int quregIds[], long numQuregs) {
    
    // check all quregs are created
    for (int i=0; i<numQuregs; i++) {
        if (!quregIsCreated[quregIds[i]]) {
            local_sendQuregNotCreatedError(quregIds[i]);
            WSPutSymbol(stdlink, "$Failed");;
            return;
        }
    }
    
    // store complex matrix as 2 flat real arrays
    long len = numQuregs * numQuregs;
    qreal* matrRe = (qreal*) malloc(len * sizeof *matrRe);
    qreal* matrIm = (qreal*) malloc(len * sizeof *matrIm);
    
    for (int r=0; r<numQuregs; r++) {
        for (int c=0; c<numQuregs; c++) {
            if (c >= r) {
                Complex val = calcInnerProduct(quregs[quregIds[r]], quregs[quregIds[c]]);
                matrRe[r*numQuregs + c] = val.real;
                matrIm[r*numQuregs + c] = val.imag;
            } else {
                matrRe[r*numQuregs + c] =   matrRe[c*numQuregs + r];
                matrIm[r*numQuregs + c] = - matrIm[c*numQuregs + r]; // conjugate transpose
            }
        }
    }
    
    WSPutFunction(stdlink, "List", 2);
    WSPutReal64List(stdlink, matrRe, len);
    WSPutReal64List(stdlink, matrIm, len);
    
    free(matrRe);
    free(matrIm);
}

/**
 * puts a Qureg into MMA, with the structure of
 * {numQubits, isDensityMatrix, realAmps, imagAmps}.
 * Instead gives -1 if error (e.g. qureg id is wrong)
 */
void internal_getStateVec(int id) {
    
    if (!local_isQuregCreated(id)) {
        local_sendQuregNotCreatedError(id);
        WSPutInteger(stdlink, -1);
        return;
    }

    Qureg qureg = quregs[id];
    syncQuESTEnv(env);
    copyStateFromGPU(qureg); // does nothing on CPU
    
    WSPutFunction(stdlink, "List", 4);
    WSPutInteger(stdlink, qureg.numQubitsRepresented);
    WSPutInteger(stdlink, qureg.isDensityMatrix);
    WSPutReal64List(stdlink, qureg.stateVec.real, qureg.numAmpsTotal);
    WSPutReal64List(stdlink, qureg.stateVec.imag, qureg.numAmpsTotal);
}

void internal_setWeightedQureg(
    double facRe1, double facIm1, int qureg1, 
    double facRe2, double facIm2, int qureg2, 
    double facReOut, double facImOut, int outID
) {
    if (!quregIsCreated[qureg1])
        local_sendQuregNotCreatedError(qureg1);
    else if (!quregIsCreated[qureg2])
        local_sendQuregNotCreatedError(qureg2);
    else if (!local_isQuregCreated(outID))
        local_sendQuregNotCreatedError(outID);
    else 
        setWeightedQureg(
            (Complex) {.real=facRe1, .imag=facIm1}, quregs[qureg1],
            (Complex) {.real=facRe2, .imag=facIm2}, quregs[qureg2],
            (Complex) {.real=facReOut, .imag=facImOut}, quregs[outID]);
    WSPutInteger(stdlink, outID);
}


/* Evaluates the expected value of a Pauli product */
void internal_calcExpecPauliProd(int quregId, int workspaceId) {
    
    // get arguments from MMA link
    int numPaulis;
    int *pauliIntCodes;
    WSGetInteger32List(stdlink, &pauliIntCodes, &numPaulis);
    int *targs;
    WSGetInteger32List(stdlink, &targs, &numPaulis);
    
    // ensure quregs exists
    Qureg qureg = quregs[quregId];
    if (!local_isQuregCreated(quregId)) {
        local_sendQuregNotCreatedError(quregId);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    Qureg workspace = quregs[workspaceId];
    if (!local_isQuregCreated(workspaceId)) {
        local_sendQuregNotCreatedError(workspaceId);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    
    // safely cast pauli codes
    enum pauliOpType pauliCodes[numPaulis];
    for (int i=0; i<numPaulis; i++)
        pauliCodes[i] = (pauliOpType) pauliIntCodes[i];
    
    WSPutReal64(stdlink, calcExpecPauliProd(qureg, targs, pauliCodes, numPaulis, workspace));
}

void local_loadPauliSumFromMMA(int numQb, int* numTerms, enum pauliOpType** arrPaulis, qreal** termCoeffs) {
    
    // get arguments from MMA link
    int numPaulis;
    WSGetReal64List(stdlink, termCoeffs, numTerms);
    int* allPauliCodes;
    WSGetInteger32List(stdlink, &allPauliCodes, &numPaulis);
    int* allPauliTargets;
    WSGetInteger32List(stdlink, &allPauliTargets, &numPaulis);
    int* numPaulisPerTerm;
    WSGetInteger32List(stdlink, &numPaulisPerTerm, numTerms);
    
    // convert {allPauliCodes}, {allPauliTargets}, {numPaulisPerTerm}, and
    // qureg.numQubitsRepresented into {pauli-code-for-every-qubit}
    int arrLen = *numTerms * numQb;
    *arrPaulis = (pauliOpType*) malloc(arrLen * sizeof **arrPaulis);
    for (int i=0; i < arrLen; i++)
        (*arrPaulis)[i] = PAULI_I;
    
    int allPaulisInd = 0;
    for (int t=0;  t < *numTerms; t++) {
        for (int j=0; j < numPaulisPerTerm[t]; j++) {
            int arrInd = t*numQb + allPauliTargets[allPaulisInd];
            (*arrPaulis)[arrInd] = (pauliOpType) allPauliCodes[allPaulisInd++];
        }
    }
    
    WSReleaseInteger32List(stdlink, allPauliCodes, numPaulis);
    WSReleaseInteger32List(stdlink, allPauliTargets, numPaulis);
    // arrPaulis and termCoeffs must be freed by caller using free and WSDisown respectively
}

void internal_calcExpecPauliSum(int quregId, int workspaceId) {
        
    // ensure quregs exists
    Qureg qureg = quregs[quregId];
    if (!local_isQuregCreated(quregId)) {
        local_sendQuregNotCreatedError(quregId);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    Qureg workspace = quregs[workspaceId];
    if (!local_isQuregCreated(workspaceId)) {
        local_sendQuregNotCreatedError(workspaceId);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    
    int numTerms; enum pauliOpType* arrPaulis; qreal* termCoeffs;
    local_loadPauliSumFromMMA(qureg.numQubitsRepresented, &numTerms, &arrPaulis, &termCoeffs);
    
    qreal val = calcExpecPauliSum(qureg, arrPaulis, termCoeffs, numTerms, workspace);
    
    free(arrPaulis);
    WSReleaseReal64List(stdlink, termCoeffs, numTerms);
    
    WSPutReal64(stdlink, val);
}

void internal_calcPauliSumMatrix(int numQubits) {
    
    int numTerms; enum pauliOpType* arrPaulis; qreal* termCoeffs;
    local_loadPauliSumFromMMA(numQubits, &numTerms, &arrPaulis, &termCoeffs);

    // create states needed to apply Pauli products
    Qureg inQureg = createQureg(numQubits, env);
    Qureg outQureg = createQureg(numQubits, env);

    // get result of paulis on each basis state
    long long int dim = inQureg.numAmpsTotal;
    WSPutFunction(stdlink, "List", 2*dim);
    
    for (long long int i=0; i < dim; i++) {
        initClassicalState(inQureg, i);
        applyPauliSum(inQureg, arrPaulis, termCoeffs, numTerms, outQureg);
        
        syncQuESTEnv(env);
        copyStateFromGPU(outQureg); // does nothing on CPU
        
        WSPutReal64List(stdlink, outQureg.stateVec.real, dim);
        WSPutReal64List(stdlink, outQureg.stateVec.imag, dim);
    }    
    
    destroyQureg(inQureg, env);
    destroyQureg(outQureg, env);
    free(arrPaulis);
    WSReleaseReal64List(stdlink, termCoeffs, numTerms);
}

void internal_applyPauliSum(int inId, int outId) {

    // ensure quregs exists
    Qureg inQureg = quregs[inId];
    if (!local_isQuregCreated(inId)) {
        local_sendQuregNotCreatedError(inId);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    Qureg outQureg = quregs[outId];
    if (!local_isQuregCreated(outId)) {
        local_sendQuregNotCreatedError(outId);
        WSPutFunction(stdlink, "Abort", 0);
        return;
    }
    
    int numTerms; enum pauliOpType* arrPaulis; qreal* termCoeffs;
    local_loadPauliSumFromMMA(inQureg.numQubitsRepresented, &numTerms, &arrPaulis, &termCoeffs);
    
    applyPauliSum(inQureg, arrPaulis, termCoeffs, numTerms, outQureg);
    
    free(arrPaulis);
    WSReleaseReal64List(stdlink, termCoeffs, numTerms);
    WSPutInteger(stdlink, outId);
}


/**
 * Frees all quregs
 */
void callable_destroyAllQuregs(void) {
    
    for (size_t id=0; id < quregs.size(); id++) {
        if (local_isQuregCreated(id)) {
            destroyQureg(quregs[id], env);
            quregIsCreated[id] = false;
        }
    }
    WSPutInteger(stdlink, -1);
}

/** 
 * Returns a list of all created quregs
 */
void callable_getAllQuregs(void) {
    
    // collect all created quregs
    int numQuregs = 0;
    int idList[quregs.size()];
    for (size_t id=0; id < quregs.size(); id++)
        if (quregIsCreated[id])
            idList[numQuregs++] = id;
    
    WSPutIntegerList(stdlink, idList, numQuregs);
}



int main(int argc, char* argv[]) {
    
    // create the single, global QuEST execution env
    env = createQuESTEnv();
    
    // establish link with MMA
	return WSMain(argc, argv);
}
