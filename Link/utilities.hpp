
#ifndef UTILITIES_H
#define UTILITIES_H

#include "QuEST.h"
#include "QuEST_complex.h"

#include <vector>



typedef std::vector<qcomp> qvector;

typedef std::vector<std::vector<qcomp>> qmatrix;

qvector local_getQvector(int dim);

qmatrix local_getQmatrix(int dim);

qmatrix operator + (const qmatrix& m1, const qmatrix& m2);


long long int local_getNumScalarsToFormMatrix(int numQubits);

qmatrix local_getQmatrixFromFlatList(qreal* flatElems, int dim);

qmatrix local_getKrausSuperoperatorFromFlatList(qreal* flatElems, int numQubits);


ComplexMatrix2 local_getMatrix2FromFlatList(qreal* flatElems);

ComplexMatrix4 local_getMatrix4FromFlatList(qreal* flatElems);

void local_setMatrixNFromFlatList(qreal* list, ComplexMatrixN m, int numQubits);

void local_createManyMatrixNFromFlatList(qreal* list, ComplexMatrixN* matrs, int numOps, int numQubits);


void local_setFlatListFromMatrixN(qreal* list, ComplexMatrixN m, int numQubits);

void local_setFlatListToMatrixDagger(qreal* list, int numQubits);

void local_setFlatListFromQmatrix(qreal* list, qmatrix m);

void local_setMatrixNFromQmatrix(ComplexMatrixN cm, qmatrix qm);


bool local_isInvertible(qmatrix matr);

bool local_isNonZero(qreal scalar);

qmatrix local_getInverse(qmatrix matr);



#endif // UTILITIES_H