#ifndef DECODERS_H
#define DECODERS_H

#include "QuEST.h"
#include "QuEST_complex.h"

#include <string>

#include "utilities.hpp"



std::string local_getCommaSep(int*   elems, int len);
std::string local_getCommaSep(qreal* elems, int len);
std::string local_getCommaSep(qcomp* elems, int len);

std::string local_qcompToStr(qcomp s);
std::string local_qmatrixToStr(qmatrix m);

std::string local_getStandardFormFromMMA(std::string expr);

void local_sendMatrixToMMA(qmatrix matrix);

void local_loadEncodedPauliStringFromMMA(
    int* numPaulis, int* numTerms, qreal** termCoeffs, int** allPauliCodes, int** allPauliTargets, int** numPaulisPerTerm);

pauliOpType* local_decodePauliString(
    int numQb, int numTerms, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm);

void local_freePauliString(
    int numPaulis, int numTerms, qreal* termCoeffs, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm, pauliOpType* arrPaulis);

PauliHamil local_loadPauliHamilForQuregFromMMA(int numQubits);

void local_freePauliHamil(PauliHamil hamil);



# endif // DECODERS_H