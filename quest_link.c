#include "wstp.h"
#include <QuEST.h>

/**
 * Codes for specifiying circuits from MMA
 */
#define OPCODE_H 0
#define OPCODE_X 1
#define OPCODE_Y 2
#define OPCODE_Z 3
#define OPCODE_Rx 4
#define OPCODE_Ry 5
#define OPCODE_Rz 6
#define OPCODE_S 7
#define OPCODE_T 8

/**
 * Max number of quregs which can simultaneously exist
 */
# define MAX_NUM_QUREGS 1000

/**
 * Global instance of QuESTEnv, created when MMA is linked.
 */
QuESTEnv env;

/**
 * Collection of instantiated Quregs
 */
Qureg quregs[MAX_NUM_QUREGS];



/*
 * INTERNAL CODE
 * Abstracts communication between C and MMA
 */

int getNextQuregID() {
    
    for (int id=0; id < MAX_NUM_QUREGS; id++)
        if (!quregs[id].isCreated)
            return id;
    
    // @TODO: raise error; out of available quregs
    return -1;
}

int getIntFromMMA(void) {
    int num;
    WSGetInteger(stdlink, &num);
    return num;
}

qreal getRealFromMMA(void) {
    qreal num;
    WSGetReal(stdlink, &num);
    return num;
}

void sendErrorToMMA(char* err_msg) {
    WSPutFunction(stdlink, "EvaluatePacket", 1);
    WSPutFunction(stdlink, "Echo", 2);
    WSPutString(stdlink, err_msg);
    WSPutString(stdlink, "Error: ");
    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
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






/*
 * MMA INTERFACE
 * Wrappers of QuEST functions, as templated in quest_templates.tm
 */


/* qureg allocation */

int wrapper_createQureg(int numQubits) {
    int id = getNextQuregID();
    quregs[id] = createQureg(numQubits, env);
    quregs[id].isCreated = 1;
    return id;
}
int wrapper_createDensityQureg(int numQubits) {
    int id = getNextQuregID();
    quregs[id] = createDensityQureg(numQubits, env);
    quregs[id].isCreated = 1;
    return id;
}
int wrapper_destroyQureg(int id) {
    if (quregs[id].isCreated) {
        destroyQureg(quregs[id], env);
        quregs[id].isCreated = 0;
    } else {
        sendErrorToMMA("qureg has not been created.");
    }
    return id;
}


/** initial states */

int wrapper_initZeroState(int id) {
    initZeroState(quregs[id]);
    return id;
}
int wrapper_initPlusState(int id) {
    initPlusState(quregs[id]);
    return id;
}
int wrapper_initClassicalState(int id, int stateInd) {
    initClassicalState(quregs[id], stateInd);
    return id;
}
int wrapper_initPureState(int quregID, int pureID) {
    initPureState(quregs[quregID], quregs[pureID]);
    return quregID;
}
int wrapper_cloneQureg(int outID, int inID) {
    cloneQureg(quregs[outID], quregs[inID]);
    return outID;
}


/** noise */

int wrapper_applyOneQubitDephaseError(int id, int qb1, qreal prob) {
    applyOneQubitDephaseError(quregs[id], qb1, prob);
    return id;
}
int wrapper_applyTwoQubitDephaseError(int id, int qb1, int qb2, qreal prob) {
    applyTwoQubitDephaseError(quregs[id], qb1, qb2, prob);
    return id;
}
int wrapper_applyOneQubitDepolariseError(int id, int qb1, qreal prob) {
    applyOneQubitDepolariseError(quregs[id], qb1, prob);
    return id;
}
int wrapper_applyTwoQubitDepolariseError(int id, int qb1, int qb2, qreal prob) {
    applyTwoQubitDepolariseError(quregs[id], qb1, qb2, prob);
    return id;
}


/* calculations */

qreal wrapper_calcProbOfOutcome(int id, int qb, int outcome) {
    return calcProbOfOutcome(quregs[id], qb, outcome);
}
qreal wrapper_calcFidelity(int id1, int id2) {
    return calcFidelity(quregs[id1], quregs[id2]);
}


/* circuit execution */

int internal_applyCircuit(int id) {
    
    // get arguments
    Qureg qureg = quregs[id];
    int numOps;
    int *opcodes;
    int *ctrls;
    int *targs;
    qreal *params;
    WSGetInteger32List(stdlink, &opcodes, &numOps);
    WSGetInteger32List(stdlink, &ctrls, &numOps);
    WSGetInteger32List(stdlink, &targs, &numOps);
    WSGetReal64List(stdlink, &params, &numOps);
    
    // backup of initial state in case of abort
    Qureg backup;
    if (qureg.isDensityMatrix)
        backup = createDensityQureg(qureg.numQubitsRepresented, env);
    else
        backup = createQureg(qureg.numQubitsRepresented, env);
    cloneQureg(backup, qureg);
    
    // apply each op
    for (int i=0; i < numOps; i++) {
        
        // check whether the user has tried to abort
        //if (WSAbort) { // why does WSAbort not exist??
        if (WSMessageReady(stdlink)) {
            int code, arg;
            WSGetMessage(stdlink, &code, &arg);
            if (code == WSTerminateMessage || 
                code == WSInterruptMessage || 
                code == WSAbortMessage ||
                code == WSImDyingMessage) {
            
                // restore the backup
                cloneQureg(qureg, backup);
                destroyQureg(backup, env);
                return id;
            }
        }
        
        // get gate info
        int op = opcodes[i];
        int ctrl = ctrls[i];
        int targ = targs[i];
        qreal param = params[i];
        
        // no avail controls
        if (op == OPCODE_H)
            hadamard(qureg, targ);
        if (op == OPCODE_S)
            sGate(qureg, targ);
        if (op == OPCODE_T)
            tGate(qureg, targ);
        
        // no avail params
        if (op == OPCODE_X) {
            if (ctrl == -1)
                pauliX(qureg, targ);
            else
                controlledNot(qureg, ctrl, targ);
        }
        if (op == OPCODE_Y) {
            if (ctrl == -1)
                pauliY(qureg, targ);
            else
                controlledPauliY(qureg, ctrl, targ);
        }
        if (op == OPCODE_Z) {
            if (ctrl == -1)
                pauliZ(qureg, targ);
            else
                controlledPhaseFlip(qureg, ctrl, targ);
        }
        
        // params and ctrl avail
        if (op == OPCODE_Rx) {
            if (ctrl == -1)
                rotateX(qureg, targ, param);
            else
                controlledRotateX(qureg, ctrl, targ, param);
        }
        if (op == OPCODE_Ry) {
            if (ctrl == -1)
                rotateY(qureg, targ, param);
            else
                controlledRotateY(qureg, ctrl, targ, param);
        }
        if (op == OPCODE_Rz) {
            if (ctrl == -1)
                rotateZ(qureg, targ, param);
            else
                controlledRotateZ(qureg, ctrl, targ, param);
        }
        
        // @TODO must warn user if encountered invalid gate
        // @TODO (not though we can refine MMA match pattern)
    }
    
    destroyQureg(backup, env);
    return id;
}

/**
 * puts a Qureg into MMA, with the structure of
 * {numQubits, isDensityMatrix, realAmps, imagAmps}
 */
void internal_getStateVec(int id) {

    Qureg qureg = quregs[id];
    syncQuESTEnv(env);
    copyStateFromGPU(qureg); // does nothing on CPU
    
    WSPutFunction(stdlink, "List", 4);
    WSPutInteger(stdlink, qureg.numQubitsRepresented);
    WSPutInteger(stdlink, qureg.isDensityMatrix);
    WSPutReal64List(stdlink, qureg.stateVec.real, qureg.numAmpsTotal);
    WSPutReal64List(stdlink, qureg.stateVec.imag, qureg.numAmpsTotal);
}



int main(int argc, char* argv[]) {
    
    // create the single, global QuEST execution env
    env = createQuESTEnv();
    
    // indicate that no quregs have yet been created
    for (int id=0; id < MAX_NUM_QUREGS; id++)
        quregs[id].isCreated = 0;
    
    // establish link with MMA
	return WSMain(argc, argv);
}
