
#ifndef ERRORS_H
#define ERRORS_H

#include "QuEST_precision.h"

#include <string>
#include <exception>


class QuESTException : public std::exception {
public:
    std::string thrower;
    std::string message;
    QuESTException(std::string func, std::string msg, ...) {
        thrower = func;
        message = msg;
    }
};

extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc);

void local_sendErrorAndWait(std::string funcName, std::string errMsg);

void local_sendErrorAndFail(std::string funcName, std::string errMsg);

void local_sendErrorAndAbort(std::string funcName, std::string errMsg);

void local_sendErrorAndFailOrAbortFromExcep(std::string funcName, std::string excepThrower, std::string errMsg);

void local_sendWarningAndContinue(std::string funcName, std::string warnMsg);

void local_throwExcepIfQuregNotCreated(int id);

void local_throwExcepIfUserAborted();

QuESTException local_gateUnsupportedExcep(std::string gateSyntax, std::string gateName);

QuESTException local_wrongNumGateParamsExcep(std::string gateSyntax, std::string gateSymb, int wrongNumParams, int rightNumParams);

QuESTException local_wrongNumGateTargsExcep(std::string gateSyntax, std::string gateSymb, int wrongNumTargs, int rightNumTargs);
QuESTException local_wrongNumGateTargsExcep(std::string gateSyntax, std::string gateSymb, int wrongNumTargs, int rightNumA, int rightNumB);

QuESTException local_wrongNumDerivParamsExcep(std::string gate, int wrongNumParams, int rightNumParams);

QuESTException local_unrecognisedGateExcep(std::string gateSyntax, int opcode, const char* caller);

QuESTException local_invalidProbExcep(std::string gateName, qreal prob, std::string maxProb);

void sendDebugEcho(std::string msg);


# endif // ERRORS_H