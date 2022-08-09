#ifndef DECODERS_H
#define DECODERS_H

#include "QuEST.h"
#include "QuEST_complex.h"

#include "utilities.hpp"



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