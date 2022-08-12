
#include "QuEST_precision.h"
#include "QuEST_complex.h"

#include "utilities.hpp"

#define MIN_NON_ZERO_EPS_FAC 1E4



/*
 * matrix formatters
 */
 
ComplexMatrix2 local_getMatrix2FromFlatList(qreal* list) {
    int dim = 2;
    ComplexMatrix2 m;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            m.real[r][c] = list[2*(dim*r+c)];
            m.imag[r][c] = list[2*(dim*r+c)+1];
        }
    return m;
}

ComplexMatrix4 local_getMatrix4FromFlatList(qreal* list) {
    int dim = 4;
    ComplexMatrix4 m;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            m.real[r][c] = list[2*(dim*r+c)];
            m.imag[r][c] = list[2*(dim*r+c)+1];
        }
    return m;
}

void local_setMatrixNFromFlatList(qreal* list, ComplexMatrixN m, int numQubits) {
    long long int dim = (1LL << numQubits);
    for (long long int r=0; r<dim; r++)
        for (long long int c=0; c<dim; c++) {
            m.real[r][c] = list[2*(dim*r+c)];
            m.imag[r][c] = list[2*(dim*r+c)+1];
        }
}

void local_setFlatListFromMatrixN(qreal* list, ComplexMatrixN m, int numQubits) {
    long long int dim = (1LL << numQubits);
    for (long long int r=0; r<dim; r++)
        for (long long int c=0; c<dim; c++) {
            list[2*(dim*r+c)]   = m.real[r][c];
            list[2*(dim*r+c)+1] = m.imag[r][c];
        }
}

void local_setFlatListFromQmatrix(qreal* list, qmatrix m) {
    size_t dim = m.size();
    for (size_t r=0; r<dim; r++)
        for (size_t c=0; c<dim; c++) {
            list[2*(dim*r+c)]   = real(m[r][c]);
            list[2*(dim*r+c)+1] = imag(m[r][c]);
        }
}

void local_setFlatListToMatrixDagger(qreal* list, int numQubits) {
    
    long long int dim = (1LL << numQubits);
    
    for (long long int r=0; r<dim; r++) {
        for (long long int c=0; c<r; c++) {
            
            qreal tmpRe = list[2*dim*r + 2*c];
            qreal tmpIm = list[2*dim*r + 2*c + 1];
            
            list[2*dim*r + 2*c]     =   list[2*dim*c + 2*r];
            list[2*dim*r + 2*c + 1] = - list[2*dim*c + 2*r + 1];
            
            list[2*dim*c + 2*r]     =   tmpRe;
            list[2*dim*c + 2*r + 1] = - tmpIm;
        }
    }
    
    for (long long int r=0; r<dim; r++)
        list[2*dim*r + 2*r + 1] *= -1;
}

void local_createManyMatrixNFromFlatList(qreal* list, ComplexMatrixN* matrs, int numOps, int numQubits) {
    long long int dim = (1LL << numQubits);
    for (int i=0; i<numOps; i++) {
        matrs[i] = createComplexMatrixN(numQubits);
        local_setMatrixNFromFlatList(&list[2*dim*dim*i], matrs[i], numQubits);
    }
}

long long int local_getNumScalarsToFormMatrix(int numQubits) {
    long long int dim = (1LL << numQubits);
    return dim*dim*2; // fac 2 for separate real and imag cmoponents
}



/* 
 * vector algebras 
 */
 
qvector local_getQvector(int dim) {
    
    // we discourage calling qvector(dim) directly, since it may mislead
    // one to believe qmatrix(dim) is valid, which will instead cause 
    // inoccuous segmentation faults
    return qvector(dim);
}



/* 
 * matrix algebras 
 */
 
 qmatrix operator + (const qmatrix& m1, const qmatrix& m2) {
     qmatrix out = m1;
     for (size_t r=0; r<m1.size(); r++)
         for (size_t c=0; c<m1.size(); c++)
             out[r][c] += m2[r][c];
     return out;
 }

qmatrix local_getQmatrix(int dim) {
    
    qmatrix matr = qmatrix(dim);
    for (int i=0; i<dim; i++)
        matr[i].resize(dim);
        
    return matr;
}

qmatrix local_getQmatrixFromFlatList(qreal* flatElems, int dim) {
    
    qmatrix matr = local_getQmatrix(dim);
        
    int n = 0;
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++) {
            matr[i][j] = qcomp(flatElems[n], flatElems[n+1]);
            n += 2;
        }
            
    return matr;
}

void local_setMatrixNFromQmatrix(ComplexMatrixN cm, qmatrix qm) {
    size_t dim = qm.size();
    for (size_t r=0; r<dim; r++)
        for (size_t c=0; c<dim; c++) {
            cm.real[r][c] = real(qm[r][c]);
            cm.imag[r][c] = imag(qm[r][c]);
        }
}

qmatrix local_getKroneckerProduct(qmatrix a, qmatrix b) {
    
    qmatrix prod = local_getQmatrix(a.size() * b.size());
    for (size_t r=0; r<b.size(); r++)
        for (size_t c=0; c<b.size(); c++)
            for (size_t i=0; i<a.size(); i++)
                for (size_t j=0; j<a.size(); j++)
                    prod[r+b.size()*i][c+b.size()*j] = a[i][j] * b[r][c];
                    
    return prod;
}

qmatrix local_getConjugate(qmatrix m) {
    
    qmatrix c = local_getQmatrix(m.size());
    for (size_t i=0; i<m.size(); i++)
        for (size_t j=0; j<m.size(); j++)
            c[i][j] = conj(m[i][j]);
            
    return c;
}

qmatrix local_getKrausSuperoperatorFromFlatList(qreal* flatElems, int numQubits) {
    
    int numOps = (int) flatElems[0];
    long long int numCols = (1LL << numQubits); // = numRows
    
    qmatrix superOp = local_getQmatrix(numCols*numCols);

    for (int n=0; n < numOps; n++) {
        qmatrix op = local_getQmatrixFromFlatList(&flatElems[1 + 2*numCols*numCols*n], numCols);
        qmatrix opConj = local_getConjugate(op);
        superOp = superOp + local_getKroneckerProduct(opConj, op);
    }
    
    return superOp;
}

qreal local_getDeterminant(qmatrix m) {
    
    if (m.size() == 2)
        return abs(m[0][0])*abs(m[1][1]) - abs(m[0][1])*abs(m[1][0]);
        
    return 0;
    
    // TODO!
}

bool local_isNonZero(qreal scalar) {
    
    return (abs(scalar) > MIN_NON_ZERO_EPS_FAC * REAL_EPS);
}

bool local_isInvertible(qmatrix matr) {
    
    qreal det = abs(local_getDeterminant(matr));
    return local_isNonZero(det);
}

qmatrix local_getInverse(qmatrix matr) {
    
    // assume matr is numerically invertible
    
    qmatrix inv = local_getQmatrix(matr.size());
    
    // TODO!
    
    return inv;
}
