
#include <random>

#include "QuEST_precision.h"
#include "QuEST_complex.h"

#include "utilities.hpp"
#include "errors.hpp"



#define MIN_NON_ZERO_EPS_FAC 1E4



/* 
 * RNG 
 */
 
std::mt19937 randGen(std::random_device{}()); // auto-seeds
std::uniform_real_distribution<qreal> randDist(0,1);
 
int local_getRandomIndex(qreal* weights, int numInds) {
    
    //weights are assumed in [0,1] and normalised to sum to 1
    
    qreal r = randDist(randGen);
    
    qreal weightSum = 0;
    
    for (int i=0; i<numInds; i++) {
        
        weightSum += weights[i];
        if (r <= weightSum)
            return i;
    }
    
    // this should never happen, unless we got exceptionally unlucky with rounding error
    return numInds-1;
}

int local_getRandomIndex(int numInds) {
    
    qreal r = randDist(randGen);
    
    qreal weight = 1/(qreal) numInds;
    qreal weightSum = 0;

    for (int i=0; i<numInds; i++) {
        
        weightSum += weight;
        if (r <= weightSum)
            return i;
    }
    
    // this should never happen, unless we got exceptionally unlucky with rounding error
    return numInds-1;
}

void local_lazyShuffle(std::vector<int> &array) {
    
    // this is a monstrously poor shuffle; merely swapping random elements a
    // fixed number of times. It should only be used when the array doesn't *really*
    // need to be shuffled, just mucked about a little for numerical reasons,
    // in runtime-performance critical code
    
    int len = (int) array.size();
    
    for (int n=0; n<len; n++) {
        
        int i = local_getRandomIndex(len);
        int j = local_getRandomIndex(len);
        
        int tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }
}



/*
 * vector formatters 
 */
 
void local_setSubDiagonalOpFromFlatList(qreal* list, SubDiagonalOp op) {
    for (long long int i=0; i<op.numElems; i++) {
        op.real[i] = list[2*i];
        op.imag[i] = list[2*i+1];
    }
}

void local_setFlatListToDiagonalMatrixDagger(qreal* list, int numQubits) {
    long long int dim = (1LL << numQubits);
    for (long long int i=0; i<dim; i++)
        list[2*i] *= -1;
}

void local_setFlatListFromQvector(qreal* list, qvector v) {
    for (size_t i=0; i<v.size(); i++) {
        list[2*i] = real(v[i]);
        list[2*i+1] = imag(v[i]);
    }
}

long long int local_getNumRealScalarsToFormDiagonalMatrix(int numQubits) {
    long long int dim = (1LL << numQubits);
    return dim*2; // fac 2 for separate real and imag cmoponents
}



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

ComplexMatrix2 local_getMatrix2FromFlatListAtIndex(qreal* list, int n) {
    int dim = 2;
    int len = 2*dim*dim;
    ComplexMatrix2 m;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            m.real[r][c] = list[len*n + 2*(dim*r+c)];
            m.imag[r][c] = list[len*n + 2*(dim*r+c)+1];
        }
    return m;
}

ComplexMatrix4 local_getMatrix4FromFlatListAtIndex(qreal* list, int n) {
    int dim = 4;
    int len = 2*dim*dim;
    ComplexMatrix4 m;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            m.real[r][c] = list[len*n + 2*(dim*r+c)];
            m.imag[r][c] = list[len*n + 2*(dim*r+c)+1];
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

void local_setMatrixNFromFlatListAtIndex(qreal* list, ComplexMatrixN m, int numQubits, int n) {
    long long int dim = (1LL << numQubits);
    long long int len = 2*dim*dim;
    for (long long int r=0; r<dim; r++)
        for (long long int c=0; c<dim; c++) {
            m.real[r][c] = list[len*n + 2*(dim*r+c)];
            m.imag[r][c] = list[len*n + 2*(dim*r+c)+1];
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

long long int local_getNumRealScalarsToFormMatrix(int numQubits) {
    long long int dim = (1LL << numQubits);
    return dim*dim*2; // fac 2 for separate real and imag cmoponents
}

void local_setComplexMatrix2RealFactor(ComplexMatrix2 *matr, qreal fac) {
    int dim = 2;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            matr->real[r][c] *= fac;
            matr->imag[r][c] *= fac;
        }
}

void local_setComplexMatrix4RealFactor(ComplexMatrix4 *matr, qreal fac) {
    int dim = 4;
    for (int r=0; r<dim; r++)
        for (int c=0; c<dim; c++) {
            matr->real[r][c] *= fac;
            matr->imag[r][c] *= fac;
        }
}

void local_setComplexMatrixToRealFactor(ComplexMatrixN matr, qreal fac) {
    long long int dim = (1LL << matr.numQubits);
    for (long long int r=0; r<dim; r++)
        for (long long int c=0; c<dim; c++) {
            matr.real[r][c] *= fac;
            matr.imag[r][c] *= fac;
        }
}

ComplexMatrix2 local_getZeroComplexMatrix2() {
    
    ComplexMatrix2 m;
    for (int i=0; i<2; i++)
        for (int j=0; j<2; j++) {
            m.real[i][j] = 0;
            m.imag[i][j] = 0;
        }
    return m;
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

qvector local_getQvectorFromFlatList(qreal* flatElems, int dim) {
    
    qvector vec = local_getQvector(dim);
        
    int n = 0;
    for (int i=0; i<dim; i++) {
        vec[i] = qcomp(flatElems[n], flatElems[n+1]);
        n += 2;
    }
            
    return vec;
}

bool local_isInvertible(qvector diagonal) {
    
    for (size_t i=0; i<diagonal.size(); i++)
        if (! local_isNonZero(diagonal[i]) )
            return false;
            
    return true;
}

qvector local_getInverse(qvector diagonal) {
    
    qvector inv = diagonal;
    
    for (size_t i=0; i<inv.size(); i++)
        inv[i] = qcomp(1,0)/diagonal[i];
        
    return inv;
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
        
    // clear
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
            matr[i][j] = 0;
        
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
    if ((long long int) dim != (1LL<<cm.numQubits))
        throw QuESTException("", "An internal error occurred. local_setMatrixNFromQmatrix() called with incompatibly-sized matrices.");
    
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

bool local_isNonZero(qreal scalar) {
    
    return (abs(scalar) > MIN_NON_ZERO_EPS_FAC * REAL_EPS);
}

bool local_isNonZero(qcomp scalar) {
    
    return local_isNonZero(abs(scalar));
}

qmatrix local_decomposeLU(qmatrix matr, std::vector<int>& pivots, int *numPivots) {
    
    size_t dim = matr.size();
    for (size_t i=0; i<dim; i++)
        pivots[i] = i;
        
    *numPivots = 0;
        
    for (size_t i=0; i<dim; i++) {
        
        // get maximum-sized elem of ith column
        qreal maxAbs = 0;
        size_t maxInd = i;
        for (size_t k=i; k<dim; k++) {
            qreal elemAbs = abs(matr[k][i]);
            if (elemAbs > maxAbs) {
                maxAbs = elemAbs;
                maxInd = k;
            }
        }
        
        // halt if singularity detected
        if (!local_isNonZero(maxAbs))
            throw QuESTException("", "Attempted LU decomposition of a non-invertible matrix.");
            
        // perform pivot if necessary
        if (maxInd != i) {
            (*numPivots)++;
            
            // pivot rows of matr (swap row i and maxInd)
            for (size_t c=0; c<dim; c++) {
                qcomp tmp = matr[i][c];
                matr[i][c] = matr[maxInd][c];
                matr[maxInd][c] = tmp;
            }
            
            // record pivot
            size_t tmpInd = pivots[i];
            pivots[i] = pivots[maxInd];
            pivots[maxInd] = tmpInd;
        }
        
        // triangularize
        for (size_t j=i+1; j<dim; j++) {
            matr[j][i] /= matr[i][i];
            
            for (size_t k=i+1; k<dim; k++)
                matr[j][k] -= matr[j][i] * matr[i][k];
        }
    }
    
    return matr;
}

qcomp local_getDeterminant(qmatrix m) {
    
    /* This function is not currently used, because the determinant is not a 
     * reliable predictor of matrix invertibility. For instance, when matrix is that 
     * describing Depol(0,3)[.5], the determinant is 10^(-80) due to the 2^6-sized 
     * sparse matrix, yet it is precisely invertible (by quasi-probability). So instead, 
     * we will duck-type, whereby local_decomposeLU() (and hence local_getInverse())
     * will throw an exception upon encountering row degeneracy (implying non-invertibility).
     */
    
    if (m.size() == 1)
        return m[0][0];
        
    if (m.size() == 2)
        return m[0][0]*m[1][1] - m[0][1]*m[1][0];
    
    int numPivots;
    std::vector<int> pivots(m.size()); // ignored
    qmatrix lu = local_decomposeLU(m, pivots, &numPivots); // throws
    
    qcomp det = 1;
    for (size_t i=0; i<lu.size(); i++)
        det *= lu[i][i];
    
    if ((numPivots - (int) m.size()) % 2)
        det *= -1;
    
    return det;
}

bool local_isInvertible(qmatrix matr) {
    
    /* see local_getDeterminant() for rant about avoiding determinant 
     * calculation
      */

    try {
        local_getInverse(matr); // throws
        return true;
    } catch (QuESTException& err) {
        return false;
    }
}

qmatrix local_getInverse(qmatrix matr) {
    
    // get LU decomposition
    int dim = (int) matr.size();
    std::vector<int> pivots(dim);
    int numPivots; // ignored
    qmatrix lu = local_decomposeLU(matr, pivots, &numPivots); // throws
    qmatrix inv = local_getQmatrix(dim);

    // populate inverse
    for (int j=0; j<dim; j++) {
        
        for (int i=0; i<dim; i++) {
            inv[i][j] = (pivots[i] == j)? 1 : 0;
            for (int k=0; k<i; k++)
                inv[i][j] -= lu[i][k] * inv[k][j];
        }
        
        for (int i=dim-1; i>=0; i--) {
            for (int k=i+1; k<dim; k++)
                inv[i][j] -= lu[i][k] * inv[k][j];
            inv[i][j] /= lu[i][i];
        }
    }
    
    return inv;
}

qmatrix local_getDagger(qmatrix matr) {
    
    qmatrix dag = local_getQmatrix(matr.size());
    
    for (size_t i=0; i<matr.size(); i++)
        for (size_t j=0; j<matr.size(); j++)
            dag[i][j] = conj(matr[j][i]);
    
    return dag;
}



/* 
 * property checkers
 */

bool local_isInt(qreal num) {
    return num == trunc(num);
}

bool local_isEncodedVector(qreal paramDim) {
    return (paramDim == 1.0);
}

bool local_isEncodedMatrix(qreal paramDim) {
    return (paramDim == 2.0);
}

bool local_isPossiblySquareMatrix(int numFlatReals) {
    
    if (numFlatReals % 2)
        return false;
        
    int dim = round(sqrt(numFlatReals/2));
    return (2*dim*dim == numFlatReals);
}


