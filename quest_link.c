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
 * Global instance of QuESTEnv, created when MMA is linked.
 */
QuESTEnv env;





/*
 * INTERNAL CODE
 * Abstracts communication between C and MMA
 */

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





/*
 * MMA INTERFACE
 * Wrappers of QuEST functions, as templated in quest_templates.tm
 */


/** initial states */

void wrapper_initZeroState(void) {
    Qureg qureg = getQuregFromMMA();
    initZeroState(qureg);
    putQuregToMMA(qureg);
    destroyQureg(qureg, env);
}
void wrapper_initPlusState(void) {
    Qureg qureg = getQuregFromMMA();
    initPlusState(qureg);
    putQuregToMMA(qureg);
    destroyQureg(qureg, env);
}
void wrapper_initClassicalState(void) {
    Qureg qureg = getQuregFromMMA();
    int stateInd = getIntFromMMA();
    initClassicalState(qureg, stateInd);
    putQuregToMMA(qureg);
    destroyQureg(qureg, env);
}
void wrapper_initPureState(void) {
    Qureg qureg = getQuregFromMMA();
    Qureg pure = getQuregFromMMA();
    initPureState(qureg, pure);
    putQuregToMMA(qureg);
    destroyQureg(qureg, env);
    destroyQureg(pure, env);
}


/** noise */

void applyNoise(int isDepol, int isDouble) {
    
    // get args, create qureg
    Qureg qureg = getQuregFromMMA();
    int qb1 = getIntFromMMA();
    int qb2 = (isDouble)? getIntFromMMA() : -1;
    qreal prob = getRealFromMMA();
    
    // appply noise
    if (isDepol) {
        if (isDouble)
            applyTwoQubitDepolariseError(qureg, qb1, qb2, prob);
        else
            applyOneQubitDepolariseError(qureg, qb1, prob);
    } else {
        if (isDouble)
            applyTwoQubitDephaseError(qureg, qb1, qb2, prob);
        else
            applyOneQubitDephaseError(qureg, qb1, prob);
    }
    
    // give qureg to MMA
    putQuregToMMA(qureg);
    
    // clean up
    destroyQureg(qureg, env);
}

void wrapper_applyOneQubitDephaseError(void) {
    applyNoise(0, 0);
}
void wrapper_applyTwoQubitDephaseError(void) {
    applyNoise(0, 1);
}
void wrapper_applyOneQubitDepolariseError(void) {
    applyNoise(1, 0);
}
void wrapper_applyTwoQubitDepolariseError(void) {
    applyNoise(1, 1);
}


/* calculations */

qreal wrapper_calcProbOfOutcome(void) {
    
    // get args, create qureg
    Qureg qureg = getQuregFromMMA();
    int qubit = getIntFromMMA();
    int outcome = getIntFromMMA();
    
    qreal prob = calcProbOfOutcome(qureg, qubit, outcome);
    
    // clean up
    destroyQureg(qureg, env);
    
    return prob;
}

qreal wrapper_calcFidelity(void) {
    
    Qureg qureg1 = getQuregFromMMA();
    Qureg qureg2 = getQuregFromMMA();
    
    qreal fid = calcFidelity(qureg1, qureg2);
    
    destroyQureg(qureg1, env);
    destroyQureg(qureg2, env);
    
    return fid;
}


/* circuit execution */

void applyCircuit(void) {
    
    // get arguments
    int numOps;
    int *opcodes;
    int *ctrls;
    int *targs;
    qreal *params;
    WSGetInteger32List(stdlink, &opcodes, &numOps);
    WSGetInteger32List(stdlink, &ctrls, &numOps);
    WSGetInteger32List(stdlink, &targs, &numOps);
    WSGetReal64List(stdlink, &params, &numOps);
    Qureg qureg = getQuregFromMMA();
    
    // apply each op
    for (int i=0; i < numOps; i++) {
        
        // check whether the user has tried to abort
        // why the eff does WSAbort not exist?
        //if (WSAbort) {
        if (WSMessageReady(stdlink)) {
            int code, param;
            WSGetMessage(stdlink, &code, &param);
            
            if (code == WSTerminateMessage || 
                code == WSInterruptMessage || 
                code == WSAbortMessage ||
                code == WSImDyingMessage) {
            
                destroyQureg(qureg, env);
                WSClearError(stdlink);
                WSNextPacket(stdlink);
                WSNewPacket(stdlink);
                WSPutSymbol(stdlink, "$Aborted");
                return;
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
    }
    
    putQuregToMMA(qureg);
    
    destroyQureg(qureg, env);
}



int main(int argc, char* argv[])
{
    env = createQuESTEnv();
    
	return WSMain(argc, argv);
}
