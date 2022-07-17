/** @file
 * Contains CPU and multithreaded specific high-performance backend facilities,
 * which can be considered as extensions to the QuEST API.
 *
 * @author Tyson Jones
 */
 
#include "QuEST.h"
#include "QuEST_cpu_internal.h"



void extension_addAdjointToSelf(Qureg qureg) {

    long long int numTasks = qureg.numAmpsPerChunk;
    int numQubits = qureg.numQubitsRepresented;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;

    long long int k, i, j, l;
    qreal tmp;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks, numQubits, vecRe,vecIm) \
    private  (k, i, j, l, tmp)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0LL; k<numTasks; k++) {

            // |k> = |j>|i>
            i = k & ((1LL << numQubits)-1);
            j = k >> numQubits;

            if (i < j) {
                // |l> = |i>|j>
                l = (i << numQubits) | j;

                tmp = vecRe[k] + vecRe[l];
                vecRe[k] = tmp;
                vecRe[l] = tmp;

                tmp = vecIm[k];
                vecIm[k] -= vecIm[l];
                vecIm[l] -= tmp;
            }
            else if (i == j) {
                vecRe[k] *= 2;
                vecIm[k] = 0;
            }
        }  
    } 
}

void extension_applyImagFactor(Qureg qureg, qreal imagFac) {
    
    long long int numTasks = qureg.numAmpsPerChunk;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
    long long int i;
    qreal tmp;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks, vecRe,vecIm, imagFac) \
    private  (i, tmp)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (i=0LL; i<numTasks; i++) {
            tmp = vecRe[i];
            vecRe[i] = - imagFac * vecIm[i];
            vecIm[i] =   imagFac * tmp;
        }  
    }
}

void extension_mixDephasingDeriv(Qureg qureg, int targetQubit, qreal probDeriv) {
    
    long long int numTasks = qureg.numAmpsPerChunk;
    int indOffset = targetQubit + qureg.numQubitsRepresented;
    qreal halfProb = probDeriv/2;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
    long long int k;
    int s;
    qreal f;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks,targetQubit,indOffset,halfProb, vecRe,vecIm) \
    private  (k, s, f)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0LL; k<numTasks; k++) {
            s = 2*(extractBit(targetQubit, k) == extractBit(indOffset, k)) - 1;
            f = halfProb*(s-1);
            vecRe[k] *= f;
            vecIm[k] *= f;
        }  
    }
    
    extension_addAdjointToSelf(qureg);
}

void extension_mixDampingDeriv(Qureg qureg, int targ, qreal prob, qreal probDeriv) {
    
    long long int numAmps = qureg.numAmpsPerChunk;
    int conjTarg = targ + qureg.numQubitsRepresented;
    qreal c1 = probDeriv/2;
    qreal c2 = - c1 / sqrt(1 - prob);
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
    long long int i, j;
    int b1, b2;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numAmps,targ,conjTarg,c1,c2, vecRe,vecIm) \
    private  (i,j,b1,b2)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (i=0LL; i<numAmps; i++) {
            b1 = extractBit(targ, i);
            b2 = extractBit(conjTarg, i);

            // TODO: refactor without branching 
            if (b2 == 0 && b1 == 0) {
                j = flipBit(flipBit(i, targ), conjTarg);
                vecRe[i] = c1 * vecRe[j];
                vecIm[i] = c1 * vecIm[j];
            } else if (b2 == 0 && b1 == 1) {
                vecRe[i] *= c2;
                vecIm[i] *= c2;
            } else if (b2 == 1 && b1 == 0) {
                vecRe[i] = 0;
                vecIm[i] = 0;
            } else {
                vecRe[i] *= - c1;
                vecIm[i] *= - c1;
            }
        }
    }

    extension_addAdjointToSelf(qureg);
}

