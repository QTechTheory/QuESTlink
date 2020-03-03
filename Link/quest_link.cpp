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
 * User-validation occurs both in the core QuEST backend, within this file, 
 * and some within the QuESTlink.m front-end. Validation problems in the backend 
 * propogate here by a QuESTException, which mentions the throwing function 
 * (a core QuEST API function) and the error message. These exceptions are caught here, 
 * and sent to the front-end as a 'Message[Func::error]' error, using Func::error 
 * string specified in quest_templates.tm or QuESTlink.m. 
 * Sometimes, code here will catch a core exception, tweak the message, and 
 * rethrow the exception to another caller in this file.
 */


#include "wstp.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <QuEST.h>
#include <string>
#include <vector>
#include <exception>

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

/*
 * Max number of target and control qubits which can be specified 
 * for an individual gate 
 */
#define MAX_NUM_TARGS_CTRLS 100

/*
 * Global instance of QuESTEnv, created when MMA is linked.
 */
QuESTEnv env;

/*
 * Collection of instantiated Quregs
 */
std::vector<Qureg> quregs;
std::vector<bool> quregIsCreated;




/* 
 * ERROR HANDLING
 */
 
class QuESTException : public std::exception {
public:
    std::string thrower;
    std::string message;
    QuESTException(std::string func, std::string msg, ...) {
        thrower = func;
        message = msg;
    }
};

/* channel core-QuEST validation errors into catchable exceptions
 */
extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
    throw QuESTException(errFunc, errMsg);
}

/* Reports an error message to MMA without closing the pipe (more output must follow).
 * funcName must have a ::error tag defined in either quest_templates.tm or 
 * QuESTlink.m
 */
void local_sendErrorAndWait(std::string funcName, std::string errMsg) {

    // send error to Mathematica
    WSPutFunction(stdlink, "EvaluatePacket", 1);

    // Message[myFunc::errormsg, err]
    WSPutFunction(stdlink, "Message", 2);

        // myFunc::errormsg = MessageName[myFunc, "errormsg"]
        WSPutFunction(stdlink, "MessageName", 2);
        WSPutSymbol(stdlink, funcName.c_str());
        WSPutString(stdlink, "error");

        WSPutString(stdlink, errMsg.c_str());

    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
    
    // a new packet is now expected; caller MUST send something else
}
void local_sendErrorAndFail(std::string funcName, std::string errMsg) {
    local_sendErrorAndWait(funcName, errMsg);
    WSPutSymbol(stdlink, "$Failed");
    
    // this closes the pipe; no further WSPut's should follow before control flow returns
}
void local_sendErrorAndAbort(std::string funcName, std::string errMsg) {
    local_sendErrorAndWait(funcName, errMsg);
    WSPutFunction(stdlink, "Abort", 0);
    
    // this closes the pipe; no further WSPut's should follow before control flow returns
}

void local_throwExcepIfQuregNotCreated(int id) {
    if (id < 0)
        throw QuESTException("", "qureg id " + std::to_string(id) + " is invalid (must be >= 0).");
    if (id >= (int) quregs.size() || !quregIsCreated[id])
        throw QuESTException("", "qureg (with id " + std::to_string(id) + ") has not been created");
}

QuESTException local_gateUnsupportedExcep(std::string gate) {
    return QuESTException("", "the gate '" + gate + "' is not supported.");    
}

QuESTException local_wrongNumGateParamsExcep(std::string gate, int wrongNumParams, int rightNumParams) {
    return QuESTException("", 
        "the gate '" + gate + "' accepts " + std::to_string(rightNumParams) + 
        " parameters, but " + std::to_string(wrongNumParams) + " were passed.");
}

QuESTException local_wrongNumGateTargsExcep(std::string gate, int wrongNumTargs, std::string rightNumTargs) {
    // rightNumTargs is a string so that it can be multiple e.g. "1 or 2"
    return QuESTException("",
        "the gate '" + gate + "' accepts " + rightNumTargs + ", but " +
         std::to_string(wrongNumTargs) + " were passed.");
}




/* 
 * QUREG MANAGEMENT
 */

size_t local_getNextQuregID(void) {
    size_t id;
    
    // check for next id
    for (id=0; id < quregs.size(); id++)
        if (!quregIsCreated[id])
            return id;
            
    // if none are available, make more space (using a blank Qureg)
    Qureg blank;
    id = quregs.size();
    quregs.push_back(blank);
    quregIsCreated.push_back(false);
    return id;
}

void wrapper_createQureg(int numQubits) {
    try { 
        size_t id = local_getNextQuregID();
        quregs[id] = createQureg(numQubits, env); // throws
        quregIsCreated[id] = true;
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CreateQureg", err.message);
    }
}

void wrapper_createDensityQureg(int numQubits) {
    try { 
        size_t id = local_getNextQuregID();
        quregs[id] = createDensityQureg(numQubits, env); // throws
        quregIsCreated[id] = true;
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CreateDensityQureg", err.message);
    }
}

void wrapper_destroyQureg(int id) {
    try { 
        local_throwExcepIfQuregNotCreated(id); // throws
        
        destroyQureg(quregs[id], env); // throws
        quregIsCreated[id] = false;
        WSPutInteger(stdlink, id);

    } catch( QuESTException& err) {
        local_sendErrorAndFail("DestroyQureg", err.message);
    }
}

void callable_destroyAllQuregs(void) {
    
    for (size_t id=0; id < quregs.size(); id++) {
        if (quregIsCreated[id]) {
            destroyQureg(quregs[id], env);
            quregIsCreated[id] = false;
        }
    }
    WSPutSymbol(stdlink, "Null");
}

void callable_createQuregs(int numQubits, int numQuregs) {
    try { 
        if (numQuregs < 0)
            throw QuESTException("", "Invalid number of quregs. Must be >= 0."); // throws
        
        int ids[numQuregs];
        for (int i=0; i < numQuregs; i++) {
            int id = local_getNextQuregID();
            ids[i] = id;
            quregs[id] = createQureg(numQubits, env); // throws (first; no cleanup needed)
            quregIsCreated[id] = true;
        }
        WSPutIntegerList(stdlink, ids, numQuregs);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CreateQuregs", err.message);
    }
}

void callable_createDensityQuregs(int numQubits, int numQuregs) {
    try {
        if (numQuregs < 0)
            throw QuESTException("", "Invalid number of quregs. Must be >= 0."); // throws
            
        int ids[numQuregs];
        for (int i=0; i < numQuregs; i++) {
            int id = local_getNextQuregID();
            ids[i] = id;
            quregs[id] = createDensityQureg(numQubits, env); // throws (first; no cleanup needed)
            quregIsCreated[id] = true;
        }
        WSPutIntegerList(stdlink, ids, numQuregs);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CreateDensityQuregs", err.message);
    }
} 




/*
 * GETTERS
 */

void internal_getAmp(int quregID) {
    
    // get args from MMA (must do this before possible early-exit)
    long int row;
    long int col;
    WSGetInteger64(stdlink, &row);
    WSGetInteger64(stdlink, &col);
    
    try { 
        local_throwExcepIfQuregNotCreated(quregID); // throws
        
        // ensure user supplied the correct number of args for the qureg type
        Qureg qureg = quregs[quregID];
        if (qureg.isDensityMatrix && col==-1)
            throw QuESTException("", "Called on a density matrix without supplying both row and column."); // throws
        if (!qureg.isDensityMatrix && col!=-1)
            throw QuESTException("", "Called on a state-vector, yet a meaningless column index was supplied."); // throws
        
        // fetch amp (can throw internal QuEST errors)
        Complex amp;
        if (qureg.isDensityMatrix)
            amp = getDensityAmp(qureg, row, col); // throws
        else
            amp = getAmp(qureg, row); // throws
        
        // return as complex number
        WSPutFunction(stdlink, "Complex", 2);
        WSPutReal64(stdlink, amp.real);
        WSPutReal64(stdlink, amp.imag);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("GetAmp", err.message);
    }
}

void callable_isDensityMatrix(int quregID) {
    try { 
        local_throwExcepIfQuregNotCreated(quregID); // throws
        Qureg qureg = quregs[quregID];
        WSPutInteger(stdlink, qureg.isDensityMatrix);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("IsDensityMatrix", err.message);
    }
}

/* puts a Qureg into MMA, with the structure of
 * {numQubits, isDensityMatrix, realAmps, imagAmps}.
 * Instead gives -1 if error (e.g. qureg id is wrong)
 */
void internal_getQuregMatrix(int id) {
    
    // superflous try-catch, but we want to grab the quregNotCreated error message
    // (and this also keeps the pattern consistent)
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
    
        Qureg qureg = quregs[id];
        syncQuESTEnv(env);       // does nothing on local
        copyStateFromGPU(qureg); // does nothing on CPU
        
        WSPutFunction(stdlink, "List", 4);
        WSPutInteger(stdlink, qureg.numQubitsRepresented);
        WSPutInteger(stdlink, qureg.isDensityMatrix);
        WSPutReal64List(stdlink, qureg.stateVec.real, qureg.numAmpsTotal);
        WSPutReal64List(stdlink, qureg.stateVec.imag, qureg.numAmpsTotal);    
        
    } catch (QuESTException& err) {
        local_sendErrorAndFail("GetQuregMatrix", err.message);
    }
}

/* Returns a list of all created quregs
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




/*
 * STATE INITIALISATION
 */

void wrapper_initZeroState(int id) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        initZeroState(quregs[id]);
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("InitZeroState", err.message);
    }
}

void wrapper_initPlusState(int id) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        initPlusState(quregs[id]);
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("InitPlusState", err.message);
    }
}

void wrapper_initClassicalState(int id, int stateInd) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        initClassicalState(quregs[id], stateInd); // throws
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("InitClassicalState", err.message);
    }
}

void wrapper_initPureState(int quregID, int pureID) {
    try {
        local_throwExcepIfQuregNotCreated(quregID); // throws
        local_throwExcepIfQuregNotCreated(pureID); // throws
        initPureState(quregs[quregID], quregs[pureID]); // throws
        WSPutInteger(stdlink, quregID);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("InitPureState", err.message);
    }
}

void wrapper_initStateFromAmps(int quregID, qreal* reals, long l1, qreal* imags, long l2) {
    try {
        local_throwExcepIfQuregNotCreated(quregID); // throws
        
        Qureg qureg = quregs[quregID];
        if (l1 != l2 || l1 != qureg.numAmpsTotal)
            throw QuESTException("", "incorrect number of amplitudes supplied. State has not been changed."); // throws
        
        initStateFromAmps(qureg, reals, imags); // possibly throws?
        WSPutInteger(stdlink, quregID);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("InitStateFromAmps", err.message);
    }
}

void wrapper_cloneQureg(int outID, int inID) {
    try {
        local_throwExcepIfQuregNotCreated(outID); // throws
        local_throwExcepIfQuregNotCreated(inID); // throws
        cloneQureg(quregs[outID], quregs[inID]); // throws
        WSPutInteger(stdlink, outID);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CloneQureg", err.message);
    }
}




/*
 * STATE MODIFICATION
 */

void wrapper_collapseToOutcome(int id, int qb, int outcome) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        collapseToOutcome(quregs[id], qb, outcome); // throws 
        WSPutInteger(stdlink, id);
    
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CollapseToOutcome", err.message);
    }
}

void internal_setWeightedQureg(
    double facRe1, double facIm1, int qureg1, 
    double facRe2, double facIm2, int qureg2, 
    double facReOut, double facImOut, int outID
) {
    try {
        local_throwExcepIfQuregNotCreated(qureg1); // throws
        local_throwExcepIfQuregNotCreated(qureg2); // throws
        local_throwExcepIfQuregNotCreated(outID); // throws
        
        setWeightedQureg(
            (Complex) {.real=facRe1, .imag=facIm1}, quregs[qureg1],
            (Complex) {.real=facRe2, .imag=facIm2}, quregs[qureg2],
            (Complex) {.real=facReOut, .imag=facImOut}, quregs[outID]); // throws
        
        WSPutInteger(stdlink, outID);
        
    } catch (QuESTException& err) {
        local_sendErrorAndFail("SetWeightedQureg", err.message);
    }
}




/*
 * DECOHERENCE
 */

void wrapper_mixDephasing(int id, int qb1, qreal prob) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        mixDephasing(quregs[id], qb1, prob); // throws
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("MixDephasing", err.message);
    }
}

void wrapper_mixTwoQubitDephasing(int id, int qb1, int qb2, qreal prob) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        mixTwoQubitDephasing(quregs[id], qb1, qb2, prob); // throws
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("MixTwoQubitDephasing", err.message);
    }
}

void wrapper_mixDepolarising(int id, int qb1, qreal prob) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        mixDepolarising(quregs[id], qb1, prob); // throws
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("MixDepolarising", err.message);
    }
}

void wrapper_mixTwoQubitDepolarising(int id, int qb1, int qb2, qreal prob) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        mixTwoQubitDepolarising(quregs[id], qb1, qb2, prob); // throws
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("MixTwoQubitDepolarising", err.message);
    }
}

void wrapper_mixDamping(int id, int qb, qreal prob) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        mixDamping(quregs[id], qb, prob); // throws
        WSPutInteger(stdlink, id);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("MixDamping", err.message);
    }
}




/*
 * CALCULATIONS
 */

void wrapper_calcProbOfOutcome(int id, int qb, int outcome) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        qreal prob = calcProbOfOutcome(quregs[id], qb, outcome); // throws
        WSPutReal64(stdlink, prob); 
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcProbOfOutCome", err.message);
    }
}

void wrapper_calcFidelity(int id1, int id2) {
    try {
        local_throwExcepIfQuregNotCreated(id1); // throws
        local_throwExcepIfQuregNotCreated(id2); // throws
        qreal fid = calcFidelity(quregs[id1], quregs[id2]); // throws
        WSPutReal64(stdlink, fid);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcFidelity", err.message);
    }
}

void wrapper_calcInnerProduct(int id1, int id2) {
    try {
        local_throwExcepIfQuregNotCreated(id1); // throws
        local_throwExcepIfQuregNotCreated(id2); // throws
        Complex res = calcInnerProduct(quregs[id1], quregs[id2]); // throws
        WSPutFunction(stdlink, "Complex", 2);
        WSPutReal64(stdlink, res.real);
        WSPutReal64(stdlink, res.imag);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcInnerProduct", err.message);
    }
}

void wrapper_calcDensityInnerProduct(int id1, int id2) {
    try {
        local_throwExcepIfQuregNotCreated(id1); // throws
        local_throwExcepIfQuregNotCreated(id2); // throws
        qreal res = calcDensityInnerProduct(quregs[id1], quregs[id2]); // throws
        WSPutReal64(stdlink, res);
        
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcDensityInnerProduct", err.message);
    }
}

void wrapper_calcPurity(int id) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        qreal pur = calcPurity(quregs[id]); // throws
        WSPutReal64(stdlink, pur);
    
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcPurity", err.message);
    }
}

void wrapper_calcTotalProb(int id) {
    try {
        local_throwExcepIfQuregNotCreated(id); // throws
        qreal prob = calcTotalProb(quregs[id]); // throws
        WSPutReal64(stdlink, prob);
    
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcTotalProb", err.message);
    }
}

void wrapper_calcHilbertSchmidtDistance(int id1, int id2) {
    try {
        local_throwExcepIfQuregNotCreated(id1); // throws
        local_throwExcepIfQuregNotCreated(id2); // throws
        qreal dist = calcHilbertSchmidtDistance(quregs[id1], quregs[id2]); // throws
        WSPutReal64(stdlink, dist);
    
    } catch( QuESTException& err) {
        local_sendErrorAndFail("CalcHilbertSchmidtDistance", err.message);
    }
}

/* returns a real vector with ith element 
 * calcDensityInnerProduct(qureg[rhoId], qureg[omegaIds[i]]) 
 */
void internal_calcDensityInnerProductsVector(int rhoId, int omegaIds[], long numOmegas) {
    
    // init to null so we can check whether it needs cleanup
    qreal* prods = NULL;
    
    try {
        local_throwExcepIfQuregNotCreated(rhoId); // throws
        for (int i=0; i<numOmegas; i++)
            local_throwExcepIfQuregNotCreated(omegaIds[i]); // throws
        
        // calculate inner products (must free this)
        prods = (qreal*) malloc(numOmegas * sizeof *prods);
        for (int i=0; i<numOmegas; i++)
            prods[i] = calcDensityInnerProduct(quregs[rhoId], quregs[omegaIds[i]]); // throws
            
        // send result to MMA
        WSPutReal64List(stdlink, prods, numOmegas);
        
        // and clean-up
        free(prods);
    
    } catch (QuESTException& err) {
        
        // may still need to clean-up
        if (prods != NULL)
            free(prods);
        
        local_sendErrorAndFail("CalcDensityInnerProducts", err.message);
    }
}

/* returns vector with ith element <qureg[braId]|qureg[ketIds[i]]> 
 */
void internal_calcInnerProductsVector(int braId, int ketIds[], long numKets) {
    
    // init to NULL so we can check if they need clean-up later
    qreal* vecRe = NULL;
    qreal* vecIm = NULL;
    
    try {
        local_throwExcepIfQuregNotCreated(braId); // throws
        for (int i=0; i<numKets; i++)
            local_throwExcepIfQuregNotCreated(ketIds[i]); // throws
            
        // calculate inner products (must free these)
        vecRe = (qreal*) malloc(numKets * sizeof *vecRe);
        vecIm = (qreal*) malloc(numKets * sizeof *vecIm);
        
        for (int i=0; i<numKets; i++) {
            Complex val = calcInnerProduct(quregs[braId], quregs[ketIds[i]]); // throws
            vecRe[i] = val.real;
            vecIm[i] = val.imag;
        }
        
        // send result to MMA
        WSPutFunction(stdlink, "List", 2);
        WSPutReal64List(stdlink, vecRe, numKets);
        WSPutReal64List(stdlink, vecIm, numKets);
        
        // clean-up
        free(vecRe);
        free(vecIm);
        
    } catch (QuESTException& err) {
            
        // may still need to clean-up
        if (vecRe != NULL)
            free(vecRe);
        if (vecIm != NULL)
            free(vecIm);
        
        local_sendErrorAndFail("CalcInnerProducts", err.message);
    }
}

/* returns real, symmetric matrix with ith jth element calcDensityInnerProduct(quregs[quregIds[i]], qureg[qurgIds[j]]) */
void internal_calcDensityInnerProductsMatrix(int quregIds[], long numQuregs) {
    
    // init to NULL so we can later check if it needs cleanup
    qreal* matr = NULL;
    
    try {
        // check all quregs are created
        for (int i=0; i<numQuregs; i++)
            local_throwExcepIfQuregNotCreated(quregIds[i]); // throws
            
        // store real matrix as `nested pointers`
        long len = numQuregs * numQuregs;
        matr = (qreal*) malloc(len * sizeof *matr);
        
        // compute matrix (exploiting matrix is symmetric)
        for (int r=0; r<numQuregs; r++) {
            for (int c=0; c<numQuregs; c++) {
                if (c >= r) {
                    qreal val = calcDensityInnerProduct(quregs[quregIds[r]], quregs[quregIds[c]]); // throws
                    matr[r*numQuregs + c] = val;
                } else {
                    matr[r*numQuregs + c] = matr[c*numQuregs + r];
                }
            }
        }
        
        // return
        WSPutReal64List(stdlink, matr, len);
        
        // clean-up
        free(matr);
        
    } catch (QuESTException& err) {
        
        // may still need clean-up
        if (matr != NULL)
            free(matr);
            
        // send error and exit
        local_sendErrorAndFail("CalcDensityInnerProducts", err.message);
    }
}

/* returns Hermitian matrix with ith jth element <qureg[i]|qureg[j]> */
void internal_calcInnerProductsMatrix(int quregIds[], long numQuregs) {
    
    // init to NULL so we can later check if it needs cleanup
    qreal* matrRe = NULL;
    qreal* matrIm = NULL;
    
    try {
        // check all quregs are created
        for (int i=0; i<numQuregs; i++)
            local_throwExcepIfQuregNotCreated(quregIds[i]); // throws
    
        // store complex matrix as 2 flat real arrays
        long len = numQuregs * numQuregs;
        matrRe = (qreal*) malloc(len * sizeof *matrRe);
        matrIm = (qreal*) malloc(len * sizeof *matrIm);
    
        // compute elems, exploiting matrix is Hermitian
        for (int r=0; r<numQuregs; r++) {
            for (int c=0; c<numQuregs; c++) {
                if (c >= r) {
                    Complex val = calcInnerProduct(quregs[quregIds[r]], quregs[quregIds[c]]); // throws
                    matrRe[r*numQuregs + c] = val.real;
                    matrIm[r*numQuregs + c] = val.imag;
                } else {
                    matrRe[r*numQuregs + c] =   matrRe[c*numQuregs + r];
                    matrIm[r*numQuregs + c] = - matrIm[c*numQuregs + r]; // conjugate transpose
                }
            }
        }
        
        // return
        WSPutFunction(stdlink, "List", 2);
        WSPutReal64List(stdlink, matrRe, len);
        WSPutReal64List(stdlink, matrIm, len);
        
        // cleanup
        free(matrRe);
        free(matrIm);
        
    } catch (QuESTException& err) {
        
        // may still need to cleanup
        if (matrRe != NULL)
            free(matrRe);
        if (matrIm != NULL)
            free(matrIm);
            
        local_sendErrorAndFail("CalcInnerProducts", err.message);
    }
}




/*
 * CIRCUIT EXECUTION 
 */

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

/* @param mesOutcomeCache may be NULL
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
    int* mesOutcomeCache,
    int* finalCtrlInd, int* finalTargInd, int* finalParamInd
    ) {
        
    int ctrlInd = 0;
    int targInd = 0;
    int paramInd = 0;
    int mesInd = 0;
    
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
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("X", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    pauliX(qureg, targs[targInd]); // throws
                else if (numCtrls == 1)
                    controlledNot(qureg, ctrls[ctrlInd], targs[targInd]); // throws
                else {
                    ComplexMatrix2 u = {
                        .real={{0,1},{1,0}},
                        .imag={{0}}};
                    multiControlledUnitary(qureg, &ctrls[ctrlInd], numCtrls, targs[targInd], u); // throws
                }
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
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("Rx", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    rotateX(qureg, targs[targInd], params[paramInd]); // throws
                else if (numCtrls == 1)
                    controlledRotateX(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]); // throws
                else
                    throw local_gateUnsupportedExcep("multi-controlled Rotate X"); // throws
                break;
                
            case OPCODE_Ry :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Ry", numParams, 1); // throws
                if (numTargs != 1)
                    throw local_wrongNumGateTargsExcep("Ry", numTargs, "1 target"); // throws
                if (numCtrls == 0)
                    rotateY(qureg, targs[targInd], params[paramInd]); // throws
                else if (numCtrls == 1)
                    controlledRotateY(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]); // throws
                else
                    throw local_gateUnsupportedExcep("multi-controlled Rotate Y"); // throws
                break;
                
            case OPCODE_Rz :
                if (numParams != 1)
                    throw local_wrongNumGateParamsExcep("Rz", numParams, 1); // throws
                if (numCtrls > 1)
                    throw local_gateUnsupportedExcep("multi-controlled Rotate Z"); // throws
                if (numCtrls == 1 && numTargs > 1)
                    throw local_gateUnsupportedExcep("multi-controlled multi-rotateZ"); // throws
                if (numTargs == 1) {
                    if (numCtrls == 0)
                        rotateZ(qureg, targs[targInd], params[paramInd]); // throws
                    if (numCtrls == 1)
                        controlledRotateZ(qureg, ctrls[ctrlInd], targs[targInd], params[paramInd]); // throws
                } else
                    multiRotateZ(qureg, &targs[targInd], numTargs, params[paramInd]); // throws
                break;
                
            case OPCODE_R: {
                if (numCtrls != 0)
                    throw local_gateUnsupportedExcep("controlled multi-rotate-Pauli"); // throws
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
                multiRotatePauli(qureg, &targs[targInd], paulis, numTargs, params[paramInd]); // throws
            }
                break;
            
            case OPCODE_U : {
                if (numTargs == 1 && numParams != 2*2*2)
                    throw QuESTException("", "single qubit U accepts only 2x2 matrices"); // throws
                if (numTargs == 2 && numParams != 4*4*2)
                    throw QuESTException("", "two qubit U accepts only 4x4 matrices"); // throws
                if (numTargs != 1 && numTargs != 2)
                    throw local_wrongNumGateTargsExcep("U", numTargs, "1 or 2 targets"); // throws
                
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
                    ComplexMatrix2 u = {
                        .real={{0,1},{1,0}},
                        .imag={{0}}};
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
                    if (mesOutcomeCache != NULL)
                        mesOutcomeCache[mesInd++] = outcomeVal;
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
                if (numParams > 1)
                    for (int q=0; q < numParams; q++)
                        collapseToOutcome(qureg, targs[targInd+q], (int) params[paramInd+q]); // throws
                else {
                    // check value isn't impossibly high
                    if (params[paramInd] >= (1LL << numTargs))
                        throw QuESTException("",
                            "P[ " + std::to_string((int) params[paramInd]) + "] was applied to " +
                            std::to_string(numTargs) + " qubits and exceeds their maximum represented " +
                            "value of " + std::to_string(1LL << numTargs) + "."); // throws
                    // work out each bit outcome and apply; right most (least significant) bit acts on right-most target
                    for (int q=0; q < numTargs; q++)
                        collapseToOutcome(qureg, targs[targInd+numTargs-q-1], (((int) params[paramInd]) >> q) & 1); // throws
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
                        std::to_string(numKrausOps) + " operators were passed to single-qubit Kraus[ops], " + 
                        "which accepts only >0 and <=" + std::to_string((numTargs==1)? 4:16) + " operators!"); // throws
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
                    Complex zero = (Complex) {.real=0, .imag=0};
                    Complex fac = (Complex) {.real=cos(params[paramInd]), .imag=sin(params[paramInd])};
                    setWeightedQureg(zero, qureg, zero, qureg, fac, qureg); // throws
                }
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
void internal_applyCircuit(int id) {
    
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
    Qureg backup = createCloneQureg(qureg, env); // must clean-up
    
    // count the total number of measurements performed in a circuit
    int totalNumMesGates = 0;
    int totalNumMeasurements = 0;
    for (int opInd=0; opInd < numOps; opInd++)
        if (opcodes[opInd] == OPCODE_M) {
            totalNumMesGates++;
            totalNumMeasurements += numTargsPerOp[opInd];
        }
        
    // prepare records of measurement outcomes
    int* mesOutcomeCache = (int*) malloc(totalNumMeasurements * sizeof(int)); // must clean-up
    int mesInd = 0;

    // attempt to apply the circuit
    try {
        // these fields are ignored
        int finalCtrlInd, finalTargInd, finalParamInd;
        
        local_applyGates(
            qureg, numOps, opcodes, ctrls, numCtrlsPerOp, 
            targs, numTargsPerOp, params, numParamsPerOp,
            mesOutcomeCache,
            &finalCtrlInd, &finalTargInd, &finalParamInd); // throws
            
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
        
        // clean-up
        free(mesOutcomeCache);
        destroyQureg(backup, env);
        local_freeCircuit(
            opcodes, ctrls, numCtrlsPerOp, targs, 
            numTargsPerOp, params, numParamsPerOp,
            numOps, totalNumCtrls, totalNumTargs, totalNumParams);
    
    } catch (QuESTException& err) {
        
        // restore backup
        cloneQureg(qureg, backup);
        
        // all objs need cleaning
        free(mesOutcomeCache);
        destroyQureg(backup, env);
        local_freeCircuit(
            opcodes, ctrls, numCtrlsPerOp, targs, 
            numTargsPerOp, params, numParamsPerOp,
            numOps, totalNumCtrls, totalNumTargs, totalNumParams);
            
        // report error, depending on type
        std::string backupNotice = " The qureg (" + std::to_string(id) + 
            ") has been restored to its prior state.";
        
        if (err.thrower == "")
            local_sendErrorAndFail("ApplyCircuit", err.message + backupNotice);
        else if (err.thrower == "Abort")
            local_sendErrorAndAbort("ApplyCircuit", err.message + backupNotice);
        else 
            local_sendErrorAndFail("ApplyCircuit", 
                "Error in " + err.thrower + ": " + err.message + backupNotice);
    }
}

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

    // don't record measurement outcomes
    int* mesOutcomes = NULL;
    
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
            mesOutcomes, &finalCtrlInd, &finalTargInd, &finalParamInd); // throws 

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
                    applyOneQubitMatrix(qureg, targs[finalTargInd], u2); // throws
                } else if (numTargs == 2) {
                    ComplexMatrix4 u4 = local_getMatrix4FromFlatList(&unitaryDerivs[unitaryDerivInd]);
                    unitaryDerivInd += 2*4*4;
                    applyTwoQubitMatrix(qureg, targs[finalTargInd], targs[finalTargInd+1], u4); // throws
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
            projectToOne(qureg, ctrls[finalCtrlInd]); // throws
        
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
            mesOutcomes, &finalCtrlInd, &finalTargInd, &finalParamInd); // throws
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
        if (quregs[initStateId].isDensityMatrix)
            throw QuESTException("", "Density matrices are not yet supported.");
        
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
            if (quregs[quregIds[i]].isDensityMatrix)
                throw QuESTException("", "Density matrices are not yet supported.");
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




/*
 * HAMILTONIAN EVALUATION
 */

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
            arrPaulis[arrInd] = (pauliOpType) allPauliCodes[allPaulisInd++];
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

/* Evaluates the expected value of a Pauli product */
void internal_calcExpecPauliProd(int quregId, int workspaceId) {
    
    // get arguments from MMA link before validation (must be freed)
    int numPaulis;
    int *pauliIntCodes;
    WSGetInteger32List(stdlink, &pauliIntCodes, &numPaulis);
    int *targs;
    WSGetInteger32List(stdlink, &targs, &numPaulis);
    
    try {
        local_throwExcepIfQuregNotCreated(quregId); // throws 
        local_throwExcepIfQuregNotCreated(workspaceId); // throws
        
        if (quregId == workspaceId)
            throw QuESTException("", "qureg and workspace must be different quregs.");
        
        Qureg qureg = quregs[quregId];
        Qureg workspace = quregs[workspaceId];
        
        // recast pauli codes
        enum pauliOpType pauliCodes[numPaulis];
        for (int i=0; i<numPaulis; i++)
            pauliCodes[i] = (pauliOpType) pauliIntCodes[i];
            
        qreal prod = calcExpecPauliProd(qureg, targs, pauliCodes, numPaulis, workspace); // throws
        WSPutReal64(stdlink, prod);
        
        // clean-up
        WSReleaseInteger32List(stdlink, pauliIntCodes, numPaulis);
        WSReleaseInteger32List(stdlink, targs, numPaulis);
        
    } catch (QuESTException& err) {
        
        // must still clean-up
        WSReleaseInteger32List(stdlink, pauliIntCodes, numPaulis);
        WSReleaseInteger32List(stdlink, targs, numPaulis);
        
        local_sendErrorAndFail("CalcExpecPauliProd", err.message);
    }
}

void internal_calcExpecPauliSum(int quregId, int workspaceId) {
    
    // must load MMA args before validation (these must all also be freed)
    int numPaulis, numTerms;
    qreal* termCoeffs;
    int *allPauliCodes, *allPauliTargets, *numPaulisPerTerm;
    local_loadEncodedPauliSumFromMMA(
        &numPaulis, &numTerms, &termCoeffs, &allPauliCodes, &allPauliTargets, &numPaulisPerTerm);
        
    // init to null in case loading fails, to indicate no-cleanup needed
    pauliOpType* arrPaulis = NULL;
        
    try {
        // ensure quregs exist
        local_throwExcepIfQuregNotCreated(quregId); // throws
        local_throwExcepIfQuregNotCreated(workspaceId); // throws
        
        if (quregId == workspaceId)
            throw QuESTException("", "qureg and workspace must be different quregs.");
        
        Qureg qureg = quregs[quregId];
        Qureg workspace = quregs[workspaceId];
        
        // reformat MMA args into QuEST Hamil format (must be later freed)
        arrPaulis = local_decodePauliSum(
            qureg.numQubitsRepresented, numTerms, allPauliCodes, allPauliTargets, numPaulisPerTerm); // throws
        
        // compute return value
        qreal val = calcExpecPauliSum(qureg, arrPaulis, termCoeffs, numTerms, workspace); // throws
        WSPutReal64(stdlink, val);
        
        // cleanup
        local_freePauliSum(numPaulis, numTerms, 
            termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm, arrPaulis);
    
    } catch( QuESTException& err) {
        
        // must still perform cleanup to avoid memory leak
        // (arrPaulis may still be NULL, even though other data is populated)
        local_freePauliSum(numPaulis, numTerms, 
            termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm, arrPaulis);

        local_sendErrorAndFail("CalcExpecPauliSum", err.message);
        return;
    }
}

void internal_calcPauliSumMatrix(int numQubits) {
    
    // must load MMA args before validation (these must all also be freed)
    int numPaulis, numTerms;
    qreal* termCoeffs;
    int *allPauliCodes, *allPauliTargets, *numPaulisPerTerm;
    local_loadEncodedPauliSumFromMMA(
        &numPaulis, &numTerms, &termCoeffs, &allPauliCodes, &allPauliTargets, &numPaulisPerTerm);

    // create states needed to apply Pauli products
    Qureg inQureg = createQureg(numQubits, env);
    Qureg outQureg = createQureg(numQubits, env);
    
    // init to null in case loading fails, to indicate no-cleanup needed
    pauliOpType* arrPaulis = NULL;
    
    try {
        // reformat MMA args into QuEST Hamil format (must be later freed)    
        arrPaulis = local_decodePauliSum(
            numQubits, numTerms, allPauliCodes, allPauliTargets, numPaulisPerTerm); // throws
        
        // check that applyPauliSum will succeed (small price to pay; can't interrupt loop!)
        applyPauliSum(inQureg, arrPaulis, termCoeffs, numTerms, outQureg);
        
    } catch( QuESTException& err) {

        // must still perform cleanup to avoid memory leak
        local_freePauliSum(numPaulis, numTerms, 
            termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm, arrPaulis);
        destroyQureg(inQureg, env);
        destroyQureg(outQureg, env);
        
        // then exit 
        local_sendErrorAndFail("CalcPauliSumMatrix", err.message);
        return;
    }
    
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
    
    // output has already been 'put'

    // clean up
    local_freePauliSum(numPaulis, numTerms, 
        termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm, arrPaulis);
    destroyQureg(inQureg, env);
    destroyQureg(outQureg, env);
}

void internal_applyPauliSum(int inId, int outId) {
    
    // must load MMA args before validation (these must all also be freed)
    int numPaulis, numTerms;
    qreal* termCoeffs;
    int *allPauliCodes, *allPauliTargets, *numPaulisPerTerm;
    local_loadEncodedPauliSumFromMMA(
        &numPaulis, &numTerms, &termCoeffs, &allPauliCodes, &allPauliTargets, &numPaulisPerTerm);

    // init to null in case loading fails, to indicate no-cleanup needed
    pauliOpType* arrPaulis = NULL;
    
    try {
        local_throwExcepIfQuregNotCreated(inId); // throws
        local_throwExcepIfQuregNotCreated(outId); // throws
        
        if (inId == outId)
            throw QuESTException("", "inQureg and outQureg must be different quregs.");
        
        Qureg inQureg = quregs[inId];
        Qureg outQureg = quregs[outId];
        
        // reformat MMA args into QuEST Hamil format (must be later freed)
        arrPaulis = local_decodePauliSum(
            inQureg.numQubitsRepresented, numTerms, allPauliCodes, allPauliTargets, numPaulisPerTerm); // throws
        
        applyPauliSum(inQureg, arrPaulis, termCoeffs, numTerms, outQureg); // throws
        
        // cleanup
        local_freePauliSum(numPaulis, numTerms, 
            termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm, arrPaulis);
            
        // and return
        WSPutInteger(stdlink, outId);
        
    } catch( QuESTException& err) {
        
        // must still clean-up (arrPaulis may still be NULL)
        local_freePauliSum(numPaulis, numTerms, 
            termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm, arrPaulis);
            
        // and report error
        local_sendErrorAndFail("ApplyPauliSum", err.message);
    }
}




/*
 * WSTP Launch
 */

int main(int argc, char* argv[]) {
    
    // create the single, global QuEST execution env
    env = createQuESTEnv();
    
    // establish link with MMA
	return WSMain(argc, argv);
}
