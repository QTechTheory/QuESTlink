/** @file
 * Contains CPU and multithreaded specific high-performance backend facilities,
 * which can be considered as extensions to the QuEST API.
 *
 * @author Tyson Jones
 */
 
#include "QuEST.h"
#include "QuEST_validation.h"
#include "QuEST_cpu_internal.h"

#include "errors.hpp"

#include <vector>
#include <algorithm>



bool extension_isHermitian(Qureg qureg) {

    validateDensityMatrQureg(qureg, "isHermitian (internal)");

    qreal tolerance = 1E4 * REAL_EPS;

    long long int dim = 1LL << qureg.numQubitsRepresented;
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;

    long long int c, r, i, j;

    // assume Hermitian until encountering a violating amplitude
    bool isHermit = true;

    // iterate |r><c| basis states (matrix stored column-wise, so successive r are adjacent amplitudes).
    // parallelise only outer loop, since if there a fewer threads than columns, the total work is negligible anyway
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (tolerance, dim,vecRe,vecIm, isHermit) \
    private  (c, r, i, j)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (c=0; c<dim; c++) {

            // abort if a thread has already determined non-Hermitivity
            if (!isHermit)
                continue;

            for (r=c; r<dim; r++) {
                if (!isHermit)
                    break;

                // determine |i> and |j> where |r><c| ~ |c>|r> = |i>,  |j> = |c><r| (dagger element)
                i = (c*dim) | r;
                j = (r*dim) | c;

                // non-Hermitiain if real( amp[i] ) != real( amp[j] )
                if (absReal(vecRe[i] - vecRe[j]) > tolerance)
                    isHermit = false;

                // non-Hermitiain if imag( amp[i] ) != - imag( amp[j] )
                if (absReal(vecIm[i] + vecIm[j]) > tolerance)
                    isHermit = false;
            }
        }
    }

    return isHermit;
}

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

void extension_applyRealFactor(Qureg qureg, qreal realFac) {
    
    long long int numTasks = qureg.numAmpsPerChunk;
    
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
    long long int i;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numTasks, vecRe,vecIm, realFac) \
    private  (i)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (i=0LL; i<numTasks; i++) {
            vecRe[i] *= realFac;
            vecIm[i] *= realFac;
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

void extension_calcExpecPauliProdsFromClassicalShadow(
    std::vector<qreal> &prodExpecVals, long numProds,
    int* sampleBases, int* sampleOutcomes, int numQb, long numSamples,
    int* pauliCodes, int* pauliTargs, int* numPaulisPerProd,
    int numBatches
) {
    // The input arrays are total size O(numSamples*numQb + numProds*numQb), 
    // permitting creation of temporary vectors with negligible memory overhead.
    // We hence convert the inputs to bitmasks for asymptotically faster shadow 
    // processing, as suggested by Balint Koczor
    std::vector<long> pauliIndOffset(numProds);
    std::vector<unsigned long long> pauliBitseqs(numProds);
    std::vector<unsigned long long> pauliTargBitseqs(numProds);
    std::vector<unsigned long long> outcomeTargBitseqs(numProds);
    std::vector<unsigned long long> baseBitseqs(numSamples);
    std::vector<unsigned long long> outcomeBitseqs(numSamples);
    
    // serially prepare starting indices of each pauli product, in time O(numProds)
    pauliIndOffset[0] = 0;
    for (int i=1; i<numProds; i++)
        pauliIndOffset[i] = pauliIndOffset[i-1] + numPaulisPerProd[i-1];
        
    // parallel validate inputs, avoiding breaking out of OpenMP loops
    volatile bool invalid = false;
    std::string errMsg = "";
    
    // parallel validate the Pauli products, in time O(numProds*prodSize) << O(numProds*numQb)
    unsigned long long k, numTotalPaulis;
    numTotalPaulis = pauliIndOffset[numProds-1] + numPaulisPerProd[numProds-1];
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (invalid,errMsg, numTotalPaulis, pauliCodes,pauliTargs, numQb) \
    private  (k)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0; k<numTotalPaulis; k++) {
            
            if (invalid)
                continue;
            
            if (pauliCodes[k] < 0 || pauliCodes[k] > 3) {
                errMsg = "A Pauli product contained an invalid Pauli operator code (" + 
                    std::to_string(pauliCodes[k]) + "). Each Pauli operator code must be an integer " + 
                    "0, 1, 2, 3, corresponding to Pauli operators Id, X, Y, Z respectively.";
                invalid = true;
            }
                    
            if (pauliTargs[k] < 0 || pauliTargs[k] >= numQb) {
                errMsg = "A Pauli product targeted an invalid qubit (" + 
                    std::to_string(pauliTargs[k]) + "). Note that the number of qubits in " +
                    "the shadow was inferred to be " + std::to_string(numQb) + ".";
                invalid = true;
            }
        }
    }
    
    if (invalid)
        throw QuESTException("", errMsg);
    
    // parallel validate the shadow sample bases and outcomes, in time O(numSamples*numQb)
    unsigned long long numTotalSampVals = numSamples * numQb;
    invalid = false;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (invalid,errMsg, numTotalSampVals, numSamples,sampleBases,sampleOutcomes, numQb) \
    private  (k)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (k=0; k<numTotalSampVals; k++) {
            
            if (invalid)
                continue;
            
            if (sampleBases[k] < 1 || sampleBases[k] > 3) {
                errMsg = "A shadow sample contained an invalid Pauli measurement basis code (" + 
                    std::to_string(sampleBases[k]) + "). Each Pauli measurement basis code must be an integer " + 
                    "1, 2, 3, corresponding to the X, Y, Z bases respectively.";
                invalid = true;
            }
                    
            if (sampleOutcomes[k] < 0 || sampleOutcomes[k] > 1) {
                errMsg = "A shadow sample contained an invalid qubit measurement outcome (" + 
                    std::to_string(sampleOutcomes[k]) + "). Measurement outcomes must be 0 or 1.";
                invalid = true;
            }
        }
    }
    
    if (invalid)
        throw QuESTException("", errMsg);
        
    // parallel encode paulis as length 2*numQb bit sequences, and their targets 
    // as bit masks, in time  O(numProds*prodSize) << O(numProds*numQb) 
    int p;
    long i, j;
    unsigned long long seqPaulis, seqPauliTargs, seqOutcomeTargs;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numProds,numPaulisPerProd,pauliIndOffset,pauliCodes,pauliTargs, \
              pauliBitseqs,pauliTargBitseqs,outcomeTargBitseqs) \
    private  (i,j,p, seqPaulis,seqPauliTargs,seqOutcomeTargs)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (i=0L; i<numProds; i++) {
            
            seqPaulis = 00ULL;      // length 2*numQb
            seqPauliTargs = 00ULL;  // length 2*numQb
            seqOutcomeTargs = 0ULL; // length numQb
            
            for (p=0; p<numPaulisPerProd[i]; p++) {
                j = pauliIndOffset[i] + p;
                seqPaulis |= (pauliCodes[j] << 2*pauliTargs[j]);
                seqPauliTargs |= (3ULL << 2*pauliTargs[j]);
                seqOutcomeTargs |= (1ULL << pauliTargs[j]);
            }
            
            pauliBitseqs[i] = seqPaulis;
            pauliTargBitseqs[i] = seqPauliTargs;
            outcomeTargBitseqs[i] = seqOutcomeTargs;
        }
    }

    // parallel encode shadow bases and outcomes into bit sequences, in time O(numSamples*numQb)
    int q, ind;
    long s, offset;
    unsigned long long seqBases, seqOuts;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numSamples,numQb,sampleBases,sampleOutcomes, baseBitseqs,outcomeBitseqs) \
    private  (s,seqBases,seqOuts,offset,q,ind)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (s=0; s<numSamples; s++) {
            
            seqBases = 00ULL;   // length 2*numQb
            seqOuts = 0ULL;     // length numQb
            offset = numQb*s;
            
            for (q=numQb-1; q>=0; q--) {
                ind = offset+q;
                seqBases = (seqBases << 2) | sampleBases[ind];
                seqOuts = (seqOuts << 1) | sampleOutcomes[ind];
            }
            
            baseBitseqs[s] = seqBases;
            outcomeBitseqs[s] = seqOuts;
        }
    }
    
    // parallel evauate the expected values of each pauli product, in time O(numProds*numSamples),
    // modifying output array prodExpecVals.
    unsigned long long targOuts;
    long numSampsPerBatch = (long) ceil(numSamples/(qreal) numBatches);
    int numProdPaulis, b, isFinalBatch, match, par; // re-using p
    long batchSize, val; // re-using s, i
    qreal fac;
    std::vector<qreal> batchVals(numBatches); // thread-private
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (numProds,numBatches,numSamples,numSampsPerBatch, pauliCodes,pauliTargs, \
              numPaulisPerProd, baseBitseqs,outcomeBitseqs,pauliBitseqs, \
              pauliTargBitseqs,outcomeTargBitseqs, prodExpecVals) \
    private  (numProdPaulis,isFinalBatch,batchSize, targOuts, p,b,i,s, par,val,match,fac) \
    firstprivate (batchVals)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (p=0; p<numProds; p++) {
            
            numProdPaulis = numPaulisPerProd[p];
            
            // divide the sample into batches, where numBatches expected << 100
            for (b=0; b<numBatches; b++) {
                
                // batch indivisibility may mean final batch has fewer samples than others
                isFinalBatch = (b == (numBatches-1));
                batchSize = isFinalBatch * (numSamples-(numBatches-1)*numSampsPerBatch) + (1-isFinalBatch) * numSampsPerBatch;
                
                val = 0;
                for (i=0; i<batchSize; i++) {
                    s = b*batchSize + i;
                    
                    // determine whether this sample matches the pauli product
                    match = (baseBitseqs[s] & pauliTargBitseqs[p]) == pauliBitseqs[p];
                    
                    // obtain the outcomes of only the targeted qubits (rest forced to 0)
                    targOuts = (outcomeBitseqs[s] & outcomeTargBitseqs[p]);

                    // bit-twiddling hack to determine parity of number of 1s in targOuts
                    // https://graphics.stanford.edu/~seander/bithacks.html#ParityMultiply
                    targOuts ^= targOuts >> 1;
                    targOuts ^= targOuts >> 2;
                    targOuts = (targOuts & 0x1111111111111111UL) * 0x1111111111111111UL;
                    par = (targOuts >> 60) & 1;
                    
                    // contribute sample if match (without branching)
                    val += match * (1-2*par);
                }
                
                batchVals[b] = val / (qreal) batchSize;
            }
            
            // choose the median of the batch values
            fac = pow(3, numProdPaulis);
            std::sort(batchVals.begin(), batchVals.end());
            if (numBatches % 2)
                prodExpecVals[p] = fac * batchVals[numBatches/2];
            else 
                prodExpecVals[p] = fac * .5 * (batchVals[numBatches/2] + batchVals[numBatches/2 + 1]);
        }
    }
}

