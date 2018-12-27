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
 */


#include "wstp.h"
#include <stdio.h>
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


/**
 * Reports an error message to MMA without aborting
 */
void local_sendErrorToMMA(char* err_msg) {
    WSPutFunction(stdlink, "EvaluatePacket", 1);
    WSPutFunction(stdlink, "Echo", 2);
    WSPutString(stdlink, err_msg);
    WSPutString(stdlink, "Error: ");
    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
}

void local_quregNotCreatedError(int id) {
    char buffer[100];
    sprintf(buffer, "qureg (with id %d) has not been created.", id);
    local_sendErrorToMMA(buffer);
}


/* qureg allocation */

int local_getNextQuregID() {
    for (int id=0; id < MAX_NUM_QUREGS; id++)
        if (!quregs[id].isCreated)
            return id;
    
    local_sendErrorToMMA("Maximum number of quregs have been allocated!");
    return -1;
}

int wrapper_createQureg(int numQubits) {
    int id = local_getNextQuregID();
    quregs[id] = createQureg(numQubits, env);
    quregs[id].isCreated = 1;
    return id;
}
int wrapper_createDensityQureg(int numQubits) {
    int id = local_getNextQuregID();
    quregs[id] = createDensityQureg(numQubits, env);
    quregs[id].isCreated = 1;
    return id;
}
int wrapper_destroyQureg(int id) {
    if (quregs[id].isCreated) {
        destroyQureg(quregs[id], env);
        quregs[id].isCreated = 0;
    } else
        local_quregNotCreatedError(id);
    return id;
}


/** initial states */

int wrapper_initZeroState(int id) {
    if (quregs[id].isCreated)
        initZeroState(quregs[id]);
    else
        local_quregNotCreatedError(id);
    return id;
}
int wrapper_initPlusState(int id) {
    if (quregs[id].isCreated)
        initPlusState(quregs[id]);
    else
        local_quregNotCreatedError(id);
    return id;
}
int wrapper_initClassicalState(int id, int stateInd) {
    if (quregs[id].isCreated)
        initClassicalState(quregs[id], stateInd);
    else
        local_quregNotCreatedError(id);
    return id;
}
int wrapper_initPureState(int quregID, int pureID) {
    if (quregs[quregID].isCreated)
        initPureState(quregs[quregID], quregs[pureID]);
    else
        local_quregNotCreatedError(quregID);
    return quregID;
}
int wrapper_initStateFromAmps(int quregID, qreal* reals, int l1, qreal* imags, int l2) {
    Qureg qureg = quregs[quregID];
    
    if (!qureg.isCreated)
        local_quregNotCreatedError(quregID);
    else if (l1 != l2 || l1 != qureg.numAmpsTotal)
        local_sendErrorToMMA("incorrect number of amplitudes supplied.");
    else
        initStateFromAmps(qureg, reals, imags);
    return quregID;
}
int wrapper_cloneQureg(int outID, int inID) {
    if (!quregs[outID].isCreated)
        local_quregNotCreatedError(outID);
    else if (!quregs[inID].isCreated)
        local_quregNotCreatedError(inID);
    else
        cloneQureg(quregs[outID], quregs[inID]);
    return outID;
}


/** noise */

int wrapper_applyOneQubitDephaseError(int id, int qb1, qreal prob) {
    if (quregs[id].isCreated)
        applyOneQubitDephaseError(quregs[id], qb1, prob);
    else
        local_quregNotCreatedError(id);
    return id;
}
int wrapper_applyTwoQubitDephaseError(int id, int qb1, int qb2, qreal prob) {
    if (quregs[id].isCreated)
        applyTwoQubitDephaseError(quregs[id], qb1, qb2, prob);
    else
        local_quregNotCreatedError(id);
    return id;
}
int wrapper_applyOneQubitDepolariseError(int id, int qb1, qreal prob) {
    if (quregs[id].isCreated)
        applyOneQubitDepolariseError(quregs[id], qb1, prob);
    else
        local_quregNotCreatedError(id);
    return id;
}
int wrapper_applyTwoQubitDepolariseError(int id, int qb1, int qb2, qreal prob) {
    if (quregs[id].isCreated)
        applyTwoQubitDepolariseError(quregs[id], qb1, qb2, prob);
    else
        local_quregNotCreatedError(id);
    return id;
}


/* calculations */

qreal wrapper_calcProbOfOutcome(int id, int qb, int outcome) {
    if (quregs[id].isCreated)
        return calcProbOfOutcome(quregs[id], qb, outcome);
    else {
        local_quregNotCreatedError(id);
        return -1;
    }
}
qreal wrapper_calcFidelity(int id1, int id2) {
    if (!quregs[id1].isCreated) {
        local_quregNotCreatedError(id1);
        return -1;
    }
    if (!quregs[id2].isCreated) {
        local_quregNotCreatedError(id2);
        return -1;
    }
    return calcFidelity(quregs[id1], quregs[id2]);
}


/* circuit execution */

void local_backupQuregThenError(char* err_msg, int id, Qureg backup) {
    local_sendErrorToMMA(err_msg);
    cloneQureg(quregs[id], backup);
    destroyQureg(backup, env);
    
    WSPutSymbol(stdlink, "$Failed");
}

void local_gateNotValidError(char* gate, int id, Qureg backup) {
    char buffer[1000];
    sprintf(buffer, 
        "the gate '%s' is not supported. "
        "Aborting circuit and restoring qureg (id %d) to its original state.", 
        gate, id);
    local_backupQuregThenError(buffer, id, backup);
}

/** 
 * Applies a given circuit to the identified qureg.
 * The circuit is expressed as lists of opcodes (identifying gates),
 * their control qubits (-1 if not controlled), their target qubits,
 * and their parameters (0 if not parameterised).
 * The original qureg of the state is restored when this function
 * is aborted by the calling MMA, or aborted due to encountering
 * an invalid gate. In this case, $Failed is returned.
 * However, a user error caught by the QuEST backend
 * (e.g. same target and control qubit) will result in the link being
 * destroyed.
 */
void internal_applyCircuit(int id) {
    
    // @TODO: should this report superfluous/ignored params?
    // @TODO: should this error when necessary params aren't passed?
    //        wait this last is undetectable - it looks like R[0] to us
    
    // get arguments from MMA link
    int numOps, *opcodes, *ctrls, *targs;
    WSGetInteger32List(stdlink, &opcodes, &numOps);
    WSGetInteger32List(stdlink, &ctrls, &numOps);
    WSGetInteger32List(stdlink, &targs, &numOps);
    qreal *params;
    WSGetReal64List(stdlink, &params, &numOps);
    
    // ensure qureg exists
    Qureg qureg = quregs[id];
    if (!qureg.isCreated) {
        local_quregNotCreatedError(id);
        WSPutSymbol(stdlink, "$Failed");
        return;
    }
    
    // backup of initial state in case of abort
    Qureg backup;
    if (qureg.isDensityMatrix)
        backup = createDensityQureg(qureg.numQubitsRepresented, env);
    else
        backup = createQureg(qureg.numQubitsRepresented, env);
    cloneQureg(backup, qureg);
    
    // attempt to apply each gate
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
                    
                return local_backupQuregThenError(
                    "Circuit simulation aborted: restoring original qureg state.",
                    id, backup);
            }
        }

        // get gate info
        int op = opcodes[i];
        int ctrl = ctrls[i];
        int targ = targs[i];
        qreal param = params[i];
        
        switch(op) {
            
            case OPCODE_H :
                if (ctrl == -1)
                    hadamard(qureg, targ);
                else
                    return local_gateNotValidError("controlled Hadamard", id, backup);
                break;
                
            case OPCODE_S :
                if (ctrl == -1)
                    sGate(qureg, targ);
                else
                    return local_gateNotValidError("controlled S-gate", id, backup);
                break;
                
            case OPCODE_T :
                if (ctrl == -1)
                    tGate(qureg, targ);
                else
                    return local_gateNotValidError("controlled T-gate", id, backup);
                break;
        
            case OPCODE_X :
                if (ctrl == -1)
                    pauliX(qureg, targ);
                else
                    controlledNot(qureg, ctrl, targ);
                break;
                
            case OPCODE_Y :
                if (ctrl == -1)
                    pauliY(qureg, targ);
                else
                    controlledPauliY(qureg, ctrl, targ);
                break;
                
            case OPCODE_Z :
                if (ctrl == -1)
                    pauliZ(qureg, targ);
                else
                    controlledPhaseFlip(qureg, ctrl, targ);
                break;
        
            case OPCODE_Rx :
                if (ctrl == -1)
                    rotateX(qureg, targ, param);
                else
                    controlledRotateX(qureg, ctrl, targ, param);
                break;
                
            case OPCODE_Ry :
                if (ctrl == -1)
                    rotateY(qureg, targ, param);
                else
                    controlledRotateY(qureg, ctrl, targ, param);
                break;
                
            case OPCODE_Rz :
                if (ctrl == -1)
                    rotateZ(qureg, targ, param);
                else
                    controlledRotateZ(qureg, ctrl, targ, param);
                break;
                
            default:
            
                return local_backupQuregThenError(
                    "circuit contained an unknown gate (opcode: %d). "
                    "Aborting the circuit and restoring original state.",
                    id, backup);
        }
    }
    
    destroyQureg(backup, env);
    WSPutInteger(stdlink, id);
}

/**
 * puts a Qureg into MMA, with the structure of
 * {numQubits, isDensityMatrix, realAmps, imagAmps}.
 * Instead gives -1 if error (e.g. qureg id is wrong)
 */
void internal_getStateVec(int id) {
    
    if (!quregs[id].isCreated) {
        local_quregNotCreatedError(id);
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

/**
 * Frees all quregs
 */
int callable_destroyAllQuregs(void) {
    
    for (int id=0; id < MAX_NUM_QUREGS; id++) {
        if (quregs[id].isCreated) {
            destroyQureg(quregs[id], env);
            quregs[id].isCreated = 0;
        }
    }
    return -1;
}

/** 
 * Returns a list of all created quregs
 */
void callable_getAllQuregs(void) {
    
    // collect all created quregs
    int numQuregs = 0;
    int idList[MAX_NUM_QUREGS];
    for (int id=0; id < MAX_NUM_QUREGS; id++)
        if (quregs[id].isCreated)
            idList[numQuregs++] = id;
    
    WSPutIntegerList(stdlink, idList, numQuregs);
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
