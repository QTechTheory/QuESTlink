/** @file
 * Contains functions for throwing errors and propogating QuEST exceptions 
 * to Mathematica errors.
 *
 * User-validation occurs both in the core QuEST backend, within this file, 
 * and some within the QuESTlink.m front-end. Validation problems in the backend 
 * propogate here by a QuESTException, which mentions the throwing function 
 * (a core QuEST API function) and the error message. These exceptions are caught here, 
 * and sent to the front-end as a 'Message[Func::error]' error, using Func::error 
 * string specified in quest_templates.tm or QuESTlink.m. 
 * Sometimes, code QuESTlink will catch a core QuEST exception, tweak the message, and 
 * rethrow the exception to another QuESTlink catcher.
 *
 * @author Tyson Jones
 */

#include "errors.hpp"
#include "link.hpp"
#include "wstp.h"

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <exception>


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

void local_sendWarningAndContinue(std::string funcName, std::string warnMsg) {
    
    // send error to Mathematica
    WSPutFunction(stdlink, "EvaluatePacket", 1);

    // Message[myFunc::errormsg, err]
    WSPutFunction(stdlink, "Message", 2);

        // myFunc::errormsg = MessageName[myFunc, "errormsg"]
        WSPutFunction(stdlink, "MessageName", 2);
        WSPutSymbol(stdlink, funcName.c_str());
        WSPutString(stdlink, "error");              // TODO: maybe this changes too??

        WSPutString(stdlink, warnMsg.c_str());

    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
    
    // a new packet is now expected; caller MUST send something else
}

void local_sendErrorAndFailOrAbortFromExcep(std::string funcName, std::string excepThrower, std::string errMsg) {

    if (excepThrower == "")
        local_sendErrorAndFail(funcName, errMsg);
    else if (excepThrower == "Abort")
        local_sendErrorAndAbort(funcName, errMsg);
    else 
        local_sendErrorAndFail(funcName, "Cannot simulate " + excepThrower + ". " + errMsg);
}

void local_throwExcepIfQuregNotCreated(int id) {
    if (id < 0)
        throw QuESTException("", "qureg id " + std::to_string(id) + " is invalid (must be >= 0).");
    if (id >= (int) quregs.size() || !quregIsCreated[id])
        throw QuESTException("", "qureg (with id " + std::to_string(id) + ") has not been created");
}

QuESTException local_gateUnsupportedExcep(std::string gateSyntax, std::string gateName) {
    return QuESTException(gateSyntax, "The implied operation \\\"" + gateName + "\\\" is not supported.");    
} 

QuESTException local_wrongNumGateParamsExcep(std::string gateSyntax, std::string gateSymb, int wrongNumParams, int rightNumParams) {
    return QuESTException(gateSyntax, 
        "Operator " + gateSymb + " accepts " + std::to_string(rightNumParams) + " parameter" +  ((rightNumParams==1)? "":"s") + 
        ", but " + std::to_string(wrongNumParams) + " " + ((wrongNumParams==1)? "was":"were") + " passed.");
}

QuESTException local_wrongNumGateTargsExcep(std::string gateSyntax, std::string gateSymb, int wrongNumTargs, int rightNumTargs) {
    return QuESTException(gateSyntax,
        "Operator " + gateSymb + " accepts " + std::to_string(rightNumTargs) + " target qubit" + ((rightNumTargs==1)? "":"s") + 
        ", but " + std::to_string(wrongNumTargs) + " " + ((wrongNumTargs==1)? "was":"were") + " passed.");
}

QuESTException local_wrongNumGateTargsExcep(std::string gateSyntax, std::string gateSymb, int wrongNumTargs, int rightNumA, int rightNumB) {
    return QuESTException(gateSyntax,
        "Operator " + gateSymb + " accepts " + std::to_string(rightNumA) + " or " +  std::to_string(rightNumB) + 
        " target qubits, but " + std::to_string(wrongNumTargs) + " " + ((wrongNumTargs==1)? "was":"were") + " passed.");
}

QuESTException local_wrongNumDerivParamsExcep(std::string gate, int wrongNumParams, int rightNumParams) {
    return QuESTException("",
        "An internal error has occurred. The QuESTlink backend expected to receive " + std::to_string(rightNumParams) + 
        " 'derivative' scalars in order to effect the derivative of gate " + gate + ", but instead " +
        "received " + std::to_string(wrongNumParams) + ".");
}

QuESTException local_unrecognisedGateExcep(std::string gateSyntax, int opcode, const char* caller) {
    throw QuESTException(gateSyntax, 
        "An internal error has occurred. The operator (opcode " + std::to_string(opcode) + 
        ") was not recognised by internal function Gate::" + std::string(caller) + "().");
}

QuESTException local_invalidProbExcep(std::string gateName, qreal prob, std::string maxProb) {
    if (prob < 0)
        return QuESTException("", "The error probability is negative and ergo invalid.");
    return QuESTException("", 
        "The error probability exceeds the maximum of " + maxProb +
        " for which " + gateName + " is maximally mixing."); // throws
}

QuESTException local_invalidProbExcep(std::string gateName, qreal probX, qreal probY, qreal probZ) {

    if (probX < 0 || probY < 0 || probZ < 0)
        return QuESTException("", "An error probability is negative and ergo invalid.");
    if (probX > 1 || probY > 1 || probZ > 1)
        return QuESTException("", "An error probability exceeded one and is ergo invalid.");
    
    if (probX + probY + probZ > 1)
        return QuESTException("", "The error probabilities together exceeded one and are ergo invalid.");

    return QuESTException("", "The error probabilities together exceeded those which induce maximal mixing.");
}

void local_throwExcepIfUserAborted() {
    
    /* Dear ancient Wolfram Gods; why does this no longer work? 
     * Why is WSAbort undefined despite appearing in the WSTP doc?
     * Why is MLAbort undefined despite appearing in wstp.h?
     * Why is LinkSnooper not reporting a WSAbortMessage when aborting?
     */
    
    if (WSMessageReady(stdlink)) {
        int code, arg;
        WSGetMessage(stdlink, &code, &arg);
        if (code == WSTerminateMessage || code == WSInterruptMessage || 
            code == WSAbortMessage     || code == WSImDyingMessage) {
                
            throw QuESTException("Abort", "Calculation aborted."); // throws
        }
    }
}

void sendDebugEcho(std::string msg) {

    WSPutFunction(stdlink, "EvaluatePacket", 1);
    WSPutFunction(stdlink, "Echo", 1);
    WSPutString(stdlink, msg.c_str());

    WSEndPacket(stdlink);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
    
    // a new packet is now expected; caller MUST send something else
}
