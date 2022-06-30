
#ifndef ERRORS_H
#define ERRORS_H

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

void local_throwExcepIfQuregNotCreated(int id);

QuESTException local_gateUnsupportedExcep(std::string gate);

QuESTException local_wrongNumGateParamsExcep(std::string gate, int wrongNumParams, int rightNumParams);

QuESTException local_wrongNumGateTargsExcep(std::string gate, int wrongNumTargs, std::string rightNumTargs);


# endif // ERRORS_H