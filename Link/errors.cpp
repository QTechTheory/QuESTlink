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