#ifndef DECODERS_H
#define DECODERS_H

#include "QuEST.h"



ComplexMatrix2 local_getMatrix2FromFlatList(qreal* list);

ComplexMatrix4 local_getMatrix4FromFlatList(qreal* list);

void local_setMatrixNFromFlatList(qreal* list, ComplexMatrixN m, int numQubits);

void local_createManyMatrixNFromFlatList(qreal* list, ComplexMatrixN* matrs, int numOps, int numQubits);

long long int local_getNumScalarsToFormMatrix(int numQubits);


void local_loadEncodedPauliSumFromMMA(
    int* numPaulis, int* numTerms, qreal** termCoeffs, int** allPauliCodes, int** allPauliTargets, int** numPaulisPerTerm);

pauliOpType* local_decodePauliSum(
    int numQb, int numTerms, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm);

void local_freePauliSum(
    int numPaulis, int numTerms, qreal* termCoeffs, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm, pauliOpType* arrPaulis);



# endif // DECODERS_H