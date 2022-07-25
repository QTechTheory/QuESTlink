#ifndef DECODERS_H
#define DECODERS_H

#include "QuEST.h"
#include "QuEST_complex.h"



ComplexMatrix2 local_getMatrix2FromFlatList(qreal* list);

ComplexMatrix4 local_getMatrix4FromFlatList(qreal* list);

void local_setMatrixNFromFlatList(qreal* list, ComplexMatrixN m, int numQubits);

void local_createManyMatrixNFromFlatList(qreal* list, ComplexMatrixN* matrs, int numOps, int numQubits);

void local_setFlatListFromMatrixN(qreal* list, ComplexMatrixN m, int numQubits);

void local_setFlatListToMatrixDagger(qreal* list, int numQubits);

long long int local_getNumScalarsToFormMatrix(int numQubits);

void local_sendMatrixToMMA(qcomp** matrix, int dim);


void local_loadEncodedPauliStringFromMMA(
    int* numPaulis, int* numTerms, qreal** termCoeffs, int** allPauliCodes, int** allPauliTargets, int** numPaulisPerTerm);

pauliOpType* local_decodePauliString(
    int numQb, int numTerms, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm);

void local_freePauliString(
    int numPaulis, int numTerms, qreal* termCoeffs, int* allPauliCodes, int* allPauliTargets, int* numPaulisPerTerm, pauliOpType* arrPaulis);

PauliHamil local_loadPauliHamilForQuregFromMMA(int numQubits);

void local_freePauliHamil(PauliHamil hamil);



# endif // DECODERS_H