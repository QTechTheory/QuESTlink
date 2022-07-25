/** @file
 * Contains CPU and multithreaded specific high-performance backend facilities,
 * which can be considered as extensions to the QuEST API.
 *
 * @author Tyson Jones
 */
 
#include "QuEST.h"
#include "QuEST_validation.h"
#include "QuEST_cpu_internal.h"



void extension_addAdjointToSelf(Qureg qureg) {
    
    validateDensityMatrQureg(qureg, "addAdjointToSelf (internal)");

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
    
    validateDensityMatrQureg(qureg, "Deph (derivative)");
    validateTarget(qureg, targetQubit, "Deph (derivative)");
    
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

void extension_mixTwoQubitDephasingDeriv(Qureg qureg, int t1, int t2, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "two-qubit Deph (derivative)");
    validateUniqueTargets(qureg, t1, t2, "two-qubit Deph (derivative)");
    
    long long int numTasks = qureg.numAmpsPerChunk;
    int t1Shift = t1 + qureg.numQubitsRepresented;
    int t2Shift = t2 + qureg.numQubitsRepresented;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;

    qreal c = (-2/3.) * probDeriv;
    qreal f;
    long long int k;
    int b1, b2;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks, t1,t2,t1Shift,t2Shift, vecRe,vecIm, c) \
    private  (k, b1,b2, f)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0LL; k<numTasks; k++) {
            b1 = (extractBit(t1, k) == extractBit(t1Shift, k));
            b2 = (extractBit(t2, k) == extractBit(t2Shift, k));
            f =  !(b1 && b2) * c;
            vecRe[k] *= f;
            vecIm[k] *= f;
        }
    }
    
    extension_addAdjointToSelf(qureg);
}

void extension_mixDepolarisingDeriv(Qureg qureg, int targ, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "Depol (derivative)");
    validateTarget(qureg, targ, "Depol (derivative)");
    
    long long int numTasks = qureg.numAmpsPerChunk;
    int shift = qureg.numQubitsRepresented;
    int targShift = targ + shift;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
    qreal c1 = (-1/3.)*probDeriv;
    qreal c2 = 2*c1;
    
    long long int k, j;
    qreal tmpRe, tmpIm;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks, targ,shift,targShift, vecRe,vecIm, c1,c2) \
    private  (k, j, tmpRe,tmpIm)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0LL; k<numTasks; k++) {
            j = flipBit(flipBit(k, targ), targShift);

            if (extractBit(targ, k) == extractBit(targShift, k)) {
                if (j >= k) {
                    tmpRe = vecRe[k];
                    tmpIm = vecIm[k];
                    
                    vecRe[k] = c1*vecRe[k] - c1*vecRe[j];
                    vecIm[k] = c1*vecIm[k] - c1*vecIm[j];

                    vecRe[j] = c1*vecRe[j] - c1*tmpRe;
                    vecIm[j] = c1*vecIm[j] - c1*tmpIm;
                }
            } else {
                vecRe[k] *= c2;
                vecIm[k] *= c2;
            }
        }  
    }
    
    extension_addAdjointToSelf(qureg);
}

void extension_mixTwoQubitDepolarisingDeriv(Qureg qureg, int t1, int t2, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "two-qubit Depol (derivative)");
    validateUniqueTargets(qureg, t1, t2, "two-qubit Depol (derivative)");
    
    long long int numTasks = qureg.numAmpsPerChunk;
    
    int numQb = qureg.numQubitsRepresented;
    int s1 = t1 + numQb;
    int s2 = t2 + numQb;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
    qreal c1 = (-8/15.)*probDeriv;
    qreal c2 = ( 2/15.)*probDeriv;
    
    long long int k, j2, j3, j4;
    int b1, b2;
    qreal re1, re2, re3, re4;
    qreal im1, im2, im3, im4;
    qreal reSum, imSum;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks, t1,t2,s1,s2, vecRe,vecIm, c1,c2) \
    private  (k, j2,j3,j4, b1,b2, re1,re2,re3,re4,im1,im2,im3,im4, reSum,imSum)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0LL; k<numTasks; k++) {
            j2 = flipBit(flipBit(k, t1), s1);
            j3 = flipBit(flipBit(k, t2), s2);
            j4 = flipBit(flipBit(j2, t2), s2);
            
            b1 = (extractBit(t1, k) == extractBit(s1, k));
            b2 = (extractBit(t2, k) == extractBit(s2, k));
            
            if (b2 && b1) {
            
                if (k < j2 && k < j3 && k < j4) {
                    re1 = vecRe[k];     im1 = vecIm[k];
                    re2 = vecRe[j2];    im2 = vecIm[j2];
                    re3 = vecRe[j3];    im3 = vecIm[j3];
                    re4 = vecRe[j4];    im4 = vecIm[j4];
                    
                    reSum = c2 * (re1 + re2 + re3 + re4);
                    imSum = c2 * (im1 + im2 + im3 + im4);
                    
                    vecRe[k] = c1*vecRe[k] + reSum;
                    vecRe[j2] = c1*vecRe[j2] + reSum;
                    vecRe[j3] = c1*vecRe[j3] + reSum;
                    vecRe[j4] = c1*vecRe[j4] + reSum;

                    vecIm[k] = c1*vecIm[k] + imSum;
                    vecIm[j2] = c1*vecIm[j2] + imSum;
                    vecIm[j3] = c1*vecIm[j3] + imSum;
                    vecIm[j4] = c1*vecIm[j4] + imSum;
                }
            } else {
                vecRe[k] *= c1;
                vecIm[k] *= c1;
            }
        }
    }
    
    extension_addAdjointToSelf(qureg);
}

void extension_mixDampingDeriv(Qureg qureg, int targ, qreal prob, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "Damp (derivative)");
    validateTarget(qureg, targ, "Damp (derivative)");
    
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
            
            // prefer branching over cache invalidation
            if (b2 == 0 && b1 == 1) {
                vecRe[i] *= c2;
                vecIm[i] *= c2;
            } else if (b2 == 1 && b1 == 0) {
                vecRe[i] = 0;
                vecIm[i] = 0;
            } else if (b1 == 1 && b2 == 1) {
                vecRe[i] *= - c1;
                vecIm[i] *= - c1;
                
                j = flipBit(flipBit(i, targ), conjTarg);
                vecRe[j] = - vecRe[i];
                vecIm[j] = - vecIm[i];            
            }
        }
    }

    extension_addAdjointToSelf(qureg);
}

