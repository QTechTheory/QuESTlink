
#ifndef UTILITIES_H
#define UTILITIES_H

#include "QuEST.h"
#include "QuEST_complex.h"

#include <vector>



int local_getRandomIndex(qreal* weights, int numInds);

int local_getRandomIndex(int numInds);

void local_lazyShuffle(std::vector<int> &array);


typedef std::vector<qcomp> qvector;

typedef std::vector<std::vector<qcomp>> qmatrix;

qvector local_getQvector(int dim);

qmatrix local_getQmatrix(int dim);

qmatrix operator + (const qmatrix& m1, const qmatrix& m2);


long long int local_getNumRealScalarsToFormMatrix(int numQubits);

long long int local_getNumRealScalarsToFormDiagonalMatrix(int numQubits);

qvector local_getQvectorFromFlatList(qreal* flatElems, int dim);

qmatrix local_getQmatrixFromFlatList(qreal* flatElems, int dim);

qmatrix local_getKrausSuperoperatorFromFlatList(qreal* flatElems, int numQubits);


ComplexMatrix2 local_getMatrix2FromFlatList(qreal* flatElems);

ComplexMatrix4 local_getMatrix4FromFlatList(qreal* flatElems);

void local_setMatrixNFromFlatList(qreal* list, ComplexMatrixN m, int numQubits);

void local_createManyMatrixNFromFlatList(qreal* list, ComplexMatrixN* matrs, int numOps, int numQubits);

ComplexMatrix2 local_getMatrix2FromFlatListAtIndex(qreal* list, int n);

ComplexMatrix4 local_getMatrix4FromFlatListAtIndex(qreal* list, int n);

ComplexMatrix2 local_getZeroComplexMatrix2();

void local_setMatrixNFromFlatListAtIndex(qreal* list, ComplexMatrixN m, int numQubits, int n);

void local_setSubDiagonalOpFromFlatList(qreal* flatElems, SubDiagonalOp op);


void local_setFlatListFromMatrixN(qreal* list, ComplexMatrixN m, int numQubits);

void local_setFlatListToMatrixDagger(qreal* list, int numQubits);

void local_setFlatListFromQmatrix(qreal* list, qmatrix m);

void local_setFlatListFromQvector(qreal* list, qvector v);

void local_setMatrixNFromQmatrix(ComplexMatrixN cm, qmatrix qm);

void local_setComplexMatrix2RealFactor(ComplexMatrix2 *m, qreal fac);

void local_setComplexMatrix4RealFactor(ComplexMatrix4 *m, qreal fac);

void local_setComplexMatrixToRealFactor(ComplexMatrixN matr, qreal fac);

void local_setFlatListToDiagonalMatrixDagger(qreal* list, int numQubits);


bool local_isInvertible(qmatrix matr);

bool local_isInvertible(qvector diag);

bool local_isNonZero(qreal scalar);

bool local_isNonZero(qcomp scalar);

qmatrix local_getInverse(qmatrix matr); // throws

qvector local_getInverse(qvector diagonal);

qmatrix local_getDagger(qmatrix matr);


bool local_isInt(qreal num);

bool local_isEncodedVector(qreal paramDim);

bool local_isEncodedMatrix(qreal paramDim);

bool local_isPossiblySquareMatrix(int numFlatReals);


#endif // UTILITIES_H