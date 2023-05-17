/** @file
 * Contains GPU specific high-performance backend facilities,
 * which can be considered as extensions to the QuEST API.
 *
 * @author Tyson Jones
 */

#include "QuEST.h"
#include "QuEST_validation.h"

#include "errors.hpp"

#include <vector>
#include <algorithm>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>



__forceinline__ __device__ int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber) {
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

__forceinline__ __device__ long long int flipBit(const long long int number, const int bitInd) {
    return (number ^ (1LL << bitInd));
}



__global__ void extension_addAdjointToSelfKernel(Qureg qureg) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;
    
    qreal* stateRe = qureg.deviceStateVec.real;
    qreal* stateIm = qureg.deviceStateVec.imag;
    
    // |k> = |j>|i>
    int numQubits = qureg.numQubitsRepresented;
    long long int k = thisTask;
    long long int i = k & ((1LL << numQubits)-1);
    long long int j = k >> numQubits;
    
    if (i < j) {
        // |l> = |i>|j>
        long long int l = (i << numQubits) | j;

        qreal tmp = stateRe[k] + stateRe[l];
        stateRe[k] = tmp;
        stateRe[l] = tmp;

        tmp = stateIm[k];
        stateIm[k] -= stateIm[l];
        stateIm[l] -= tmp;
    }
    else if (i == j) {
        stateRe[k] *= 2;
        stateIm[k] = 0;
    }
}

void extension_addAdjointToSelf(Qureg qureg) {
    
    validateDensityMatrQureg(qureg, "addAdjointToSelf (internal)");
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_addAdjointToSelfKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg);
}



__global__ void extension_applyImagFactorKernel(Qureg qureg, qreal imagFac) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;
    
    qreal* stateRe = qureg.deviceStateVec.real;
    qreal* stateIm = qureg.deviceStateVec.imag;

    // (a + b i) (fac i) = (- fac b + a fac i)
    qreal a = stateRe[thisTask];
    qreal b = stateIm[thisTask];
    stateRe[thisTask] = - imagFac * b;
    stateIm[thisTask] = + imagFac * a;
}

void extension_applyImagFactor(Qureg qureg, qreal imagFac) {
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_applyImagFactorKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, imagFac);
}



__global__ void extension_applyRealFactorKernel(Qureg qureg, qreal realFac) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;

    qureg.deviceStateVec.real[thisTask] *= realFac;
    qureg.deviceStateVec.imag[thisTask] *= realFac;
}

void extension_applyRealFactor(Qureg qureg, qreal realFac) {
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_applyRealFactorKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, realFac);
}



__global__ void extension_mixDephasingDerivKernel(Qureg qureg, int targ, qreal probDeriv) {

    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;

    long long int targShift = targ + qureg.numQubitsRepresented;
    int s = 2*(extractBit(targ, thisTask) == extractBit(targShift, thisTask)) - 1;
    qreal f = (probDeriv/2.)*(s-1);
    qureg.deviceStateVec.real[thisTask] *= f;
    qureg.deviceStateVec.imag[thisTask] *= f;
}

void extension_mixDephasingDeriv(Qureg qureg, int targ, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "Deph (derivative)");
    validateTarget(qureg, targ, "Deph (derivative)");
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_mixDephasingDerivKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targ, probDeriv);
    extension_addAdjointToSelfKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg);
}



__global__ void extension_mixTwoQubitDephasingDerivKernel(Qureg qureg, int t1, int t2, qreal probDeriv) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;
    
    int t1Shift = t1 + qureg.numQubitsRepresented;
    int t2Shift = t2 + qureg.numQubitsRepresented;
    qreal c = (-2/3.) * probDeriv;
    
    long long int k = thisTask;
    int b1 = (extractBit(t1, k) == extractBit(t1Shift, k));
    int b2 = (extractBit(t2, k) == extractBit(t2Shift, k));
    qreal f =  !(b1 && b2) * c;
    qureg.deviceStateVec.real[k] *= f;
    qureg.deviceStateVec.imag[k] *= f;
}

void extension_mixTwoQubitDephasingDeriv(Qureg qureg, int t1, int t2, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "two-qubit Deph (derivative)");
    validateUniqueTargets(qureg, t1, t2, "two-qubit Deph (derivative)");
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_mixTwoQubitDephasingDerivKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, t1, t2, probDeriv);
    extension_addAdjointToSelfKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg);
}



__global__ void extension_mixDepolarisingDerivKernel(Qureg qureg, int targ, qreal probDeriv) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;
    
    qreal* vecRe = qureg.deviceStateVec.real;
    qreal* vecIm = qureg.deviceStateVec.imag;
    
    long long int targShift = targ + qureg.numQubitsRepresented;
    qreal c1 = (-1/3.)*probDeriv;
    qreal c2 = 2*c1;
    
    long long int k = thisTask;
    long long int j = flipBit(flipBit(k, targ), targShift);

    if (extractBit(targ, k) == extractBit(targShift, k)) {
        if (j >= k) {
            qreal tmpRe = vecRe[k];
            qreal tmpIm = vecIm[k];
            
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

void extension_mixDepolarisingDeriv(Qureg qureg, int targ, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "Depol (derivative)");
    validateTarget(qureg, targ, "Depol (derivative)");
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_mixDepolarisingDerivKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targ, probDeriv);
    extension_addAdjointToSelfKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg);
}



__global__ void extension_mixTwoQubitDepolarisingDerivKernel(Qureg qureg, int t1, int t2, qreal c1, qreal c2) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;
    
    int numQb = qureg.numQubitsRepresented;
    int s1 = t1 + numQb;
    int s2 = t2 + numQb;
    
    qreal* vecRe = qureg.deviceStateVec.real;
    qreal* vecIm = qureg.deviceStateVec.imag;
    
    long long int k = thisTask;
    long long int j2 = flipBit(flipBit(k, t1), s1);
    long long int j3 = flipBit(flipBit(k, t2), s2);
    long long int j4 = flipBit(flipBit(j2, t2), s2);
    
    int b1 = (extractBit(t1, k) == extractBit(s1, k));
    int b2 = (extractBit(t2, k) == extractBit(s2, k));
    
    if (b2 && b1) {
        if (k < j2 && k < j3 && k < j4) {
            qreal re1 = vecRe[k];     qreal im1 = vecIm[k];
            qreal re2 = vecRe[j2];    qreal im2 = vecIm[j2];
            qreal re3 = vecRe[j3];    qreal im3 = vecIm[j3];
            qreal re4 = vecRe[j4];    qreal im4 = vecIm[j4];
            
            qreal reSum = c2 * (re1 + re2 + re3 + re4);
            qreal imSum = c2 * (im1 + im2 + im3 + im4);
            
            vecRe[k]  = c1*vecRe[k]  + reSum;
            vecRe[j2] = c1*vecRe[j2] + reSum;
            vecRe[j3] = c1*vecRe[j3] + reSum;
            vecRe[j4] = c1*vecRe[j4] + reSum;

            vecIm[k]  = c1*vecIm[k]  + imSum;
            vecIm[j2] = c1*vecIm[j2] + imSum;
            vecIm[j3] = c1*vecIm[j3] + imSum;
            vecIm[j4] = c1*vecIm[j4] + imSum;
        }
    } else {
        vecRe[k] *= c1;
        vecIm[k] *= c1;
    }
}

void extension_mixTwoQubitDepolarisingDeriv(Qureg qureg, int t1, int t2, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "two-qubit Depol (derivative)");
    validateUniqueTargets(qureg, t1, t2, "two-qubit Depol (derivative)");
    
    qreal c1 = (-8/15.)*probDeriv;
    qreal c2 = ( 2/15.)*probDeriv;
    
    // 12 of every 16 amplitudes merely scale,
    // the remaining 4 are mixed
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_mixTwoQubitDepolarisingDerivKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, t1, t2, c1, c2);
    extension_addAdjointToSelfKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg);
}



__global__ void extension_mixDampingDerivKernel(Qureg qureg, int targ, qreal c1, qreal c2) {
    
    // each thread modifies one value (blegh)
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;
    
    int conjTarg = targ + qureg.numQubitsRepresented;
    qreal* vecRe = qureg.deviceStateVec.real;
    qreal* vecIm = qureg.deviceStateVec.imag;
    
    long long int i = thisTask;
    int b1 = extractBit(targ, i);
    int b2 = extractBit(conjTarg, i);

    /* easy to refactor without warp divergence, albeit with the same 
     * (possibly even worse) memory bottleneck due to non-local memory 
     * modification in the last condition 
     */
    if (b2 == 0 && b1 == 1) {
        vecRe[i] *= c2;
        vecIm[i] *= c2;
    } else if (b2 == 1 && b1 == 0) {
        vecRe[i] = 0;
        vecIm[i] = 0;
    } else if (b1 == 1 && b2 == 1){
        vecRe[i] *= - c1;
        vecIm[i] *= - c1;
        long long int j = flipBit(flipBit(i, targ), conjTarg);
        vecRe[j] = - vecRe[i];
        vecIm[j] = - vecIm[i];
    }
}

void extension_mixDampingDeriv(Qureg qureg, int targ, qreal prob, qreal probDeriv) {
    
    validateDensityMatrQureg(qureg, "Damp (derivative)");
    validateTarget(qureg, targ, "Damp (derivative)");
    
    qreal c1 = probDeriv/2.;
    qreal c2 = - c1 / sqrt(1. - prob);
    
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk/ (qreal) threadsPerCUDABlock);
    extension_mixDampingDerivKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targ, c1, c2);
    extension_addAdjointToSelfKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg);
}

__global__ void local_validateShadowPaulisKernel(int* invalid, int numQb, unsigned long long  numTotalPaulis, int* pauliCodes, int* pauliTargs) {
    long long int k = blockIdx.x*blockDim.x + threadIdx.x;
    if (k >= numTotalPaulis) return;
    
    if (pauliCodes[k] < 0 || pauliCodes[k] > 3 || pauliTargs[k] < 0 || pauliTargs[k] >= numQb)
        *invalid = 1;
}

__global__ void local_validateShadowSampsKernel(int* invalid, unsigned long long numTotalSampVals, int* sampleBases, int* sampleOutcomes) {
    long long int k = blockIdx.x*blockDim.x + threadIdx.x;
    if (k >= numTotalSampVals) return;
    
    if (sampleBases[k] < 1 || sampleBases[k] > 3 || sampleOutcomes[k] < 0 || sampleOutcomes[k] > 1)
        *invalid = 1;
}

__global__ void local_prepareShadowPauliBitseqsKernel(
    unsigned long long* pauliBitseqs, unsigned long long* pauliTargBitseqs, unsigned long long* outcomeTargBitseqs,
    long numProds, int* numPaulisPerProd, long* pauliIndOffset, int* pauliCodes, int* pauliTargs
) {
    long long int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= numProds) return;
    
    int offset = pauliIndOffset[i];
    int numPaulisInProd = numPaulisPerProd[i];
    
    unsigned long long seqPaulis = 00ULL;      // length 2*numQb
    unsigned long long seqPauliTargs = 00ULL;  // length 2*numQb
    unsigned long long seqOutcomeTargs = 0ULL; // length numQb
    
    for (int p=0; p<numPaulisInProd; p++) {
        int j = offset + p;
        seqPaulis |= (pauliCodes[j] << 2*pauliTargs[j]);
        seqPauliTargs |= (3ULL << 2*pauliTargs[j]);
        seqOutcomeTargs |= (1ULL << pauliTargs[j]);
    }
    
    pauliBitseqs[i] = seqPaulis;
    pauliTargBitseqs[i] = seqPauliTargs;
    outcomeTargBitseqs[i] = seqOutcomeTargs;
}

__global__ void local_prepareShadowSampleBitseqsKernel(
    unsigned long long* baseBitseqs, unsigned long long* outcomeBitseqs,
    int numQb, long numSamples, int* sampleBases, int* sampleOutcomes
) {
    long long int s = blockIdx.x*blockDim.x + threadIdx.x;
    if (s >= numSamples) return;
    
    long offset = numQb*s;
    unsigned long long seqBases = 00ULL;   // length 2*numQb
    unsigned long long seqOuts = 0ULL;     // length numQb

    for (int q=numQb-1; q>=0; q--) {
        int ind = offset + q;
        seqBases = (seqBases << 2) | sampleBases[ind];
        seqOuts = (seqOuts << 1) | sampleOutcomes[ind];
    }
        
    baseBitseqs[s] = seqBases;
    outcomeBitseqs[s] = seqOuts;
}

#define MAX_NUM_SHADOW_BATCHES 200
    
__global__ void extension_calcExpecPauliProdsFromClassicalShadowKernel(
    qreal* prodExpecVals,
    int numProds, int numSamples, int numBatches,
    int* numPaulisPerProd,
    unsigned long long *baseBitseqs, unsigned long long *pauliTargBitseqs, 
    unsigned long long *pauliBitseqs, unsigned long long *outcomeBitseqs, 
    unsigned long long *outcomeTargBitseqs
) {
    // parallel evauate the expected values of each pauli product, modifying output array prodExpecVals.
    long long int p = blockIdx.x*blockDim.x + threadIdx.x;
    if (p >= numProds) return;
    
    // divide this thread's job into batches 
    qreal batchVals[MAX_NUM_SHADOW_BATCHES];
    long numSampsPerBatch = (long) ceil(numSamples/(qreal) numBatches);
    int numProdPaulis = numPaulisPerProd[p];
    
    // populate batch values
    for (int b=0; b<numBatches; b++) {
        
        // batch indivisibility may mean final batch has fewer samples than others
        bool isFinalBatch = (b == (numBatches-1));
        long batchSize = isFinalBatch * (numSamples-(numBatches-1)*numSampsPerBatch) + (1-isFinalBatch) * numSampsPerBatch;
        
        long val = 0;
        for (long i=0; i<batchSize; i++) {
            long s = b*batchSize + i;
            
            // determine whether this sample matches the pauli product
            int match = (baseBitseqs[s] & pauliTargBitseqs[p]) == pauliBitseqs[p];
            
            // obtain the outcomes of only the targeted qubits (rest forced to 0)
            unsigned long long targOuts = (outcomeBitseqs[s] & outcomeTargBitseqs[p]);

            // bit-twiddling hack to determine parity of number of 1s in targOuts
            // https://graphics.stanford.edu/~seander/bithacks.html#ParityMultiply
            targOuts ^= targOuts >> 1;
            targOuts ^= targOuts >> 2;
            targOuts = (targOuts & 0x1111111111111111UL) * 0x1111111111111111UL;
            int par = (targOuts >> 60) & 1;
            
            // contribute sample if match (without branching)
            val += match * (1-2*par);
        }
        
        batchVals[b] = val / (qreal) batchSize;
    }
    
    qreal fac = pow((qreal) 3., (qreal) numProdPaulis); 
    
    // choose the median of the batch values
    thrust::sort(thrust::seq, batchVals, &(batchVals[numBatches]));
    if (numBatches % 2)
        prodExpecVals[p] = fac * batchVals[numBatches/2];
    else 
        prodExpecVals[p] = fac * .5 * (batchVals[numBatches/2] + batchVals[numBatches/2 + 1]);
}

void extension_calcExpecPauliProdsFromClassicalShadow(
    std::vector<qreal> &prodExpecVals, long numProds,
    int* sampleBases, int* sampleOutcomes, int numQb, long numSamples,
    int* pauliCodes, int* pauliTargs, int* numPaulisPerProd,
    int numBatches
) {
    // early numBatches validation, since we must be more strict than the CPU variant
    if (numBatches > MAX_NUM_SHADOW_BATCHES)
        throw QuESTException("", "The maximum number of batches permitted in GPU mode is " + std::to_string(MAX_NUM_SHADOW_BATCHES));
    
    // serial array encode: O(numProds)
    std::vector<long> pauliIndOffset(numProds);
    pauliIndOffset[0] = 0;
    for (int i=1; i<numProds; i++)
        pauliIndOffset[i] = pauliIndOffset[i-1] + numPaulisPerProd[i-1];
        
    // prepare device copy of output structure
    qreal* d_prodExpecVals;
    size_t memExpecVals = numProds * sizeof *d_prodExpecVals;
    cudaMalloc(&d_prodExpecVals, memExpecVals);
        
    // prepare device copies of input structures
    unsigned long long numTotalPaulis = pauliIndOffset[numProds-1] + numPaulisPerProd[numProds-1];
    size_t memSamps = numSamples*numQb * sizeof(int);
    size_t memPaulis = numTotalPaulis * sizeof (int);
    size_t memProds = numProds * sizeof (int);
    int* d_sampleBases;         cudaMalloc(&d_sampleBases, memSamps);       cudaMemcpy(d_sampleBases, sampleBases, memSamps, cudaMemcpyHostToDevice);
    int* d_sampleOutcomes;      cudaMalloc(&d_sampleOutcomes, memSamps);    cudaMemcpy(d_sampleOutcomes, sampleOutcomes, memSamps, cudaMemcpyHostToDevice);
    int* d_pauliCodes;          cudaMalloc(&d_pauliCodes, memPaulis);       cudaMemcpy(d_pauliCodes, pauliCodes, memPaulis, cudaMemcpyHostToDevice);
    int* d_pauliTargs;          cudaMalloc(&d_pauliTargs, memPaulis);       cudaMemcpy(d_pauliTargs, pauliTargs, memPaulis, cudaMemcpyHostToDevice);
    int* d_numPaulisPerProd;    cudaMalloc(&d_numPaulisPerProd, memProds);  cudaMemcpy(d_numPaulisPerProd, numPaulisPerProd, memProds, cudaMemcpyHostToDevice);
    
    // prepare device working structures
    size_t memOffset = numProds * sizeof(long);
    long *d_pauliIndOffset; 
    cudaMalloc(&d_pauliIndOffset, memOffset);
    cudaMemcpy(d_pauliIndOffset, pauliIndOffset.data(), memOffset, cudaMemcpyHostToDevice);
    
    memProds = numProds * sizeof(unsigned long long);
    memSamps = numSamples * sizeof(unsigned long long);
    unsigned long long *d_pauliBitseqs;         cudaMalloc(&d_pauliBitseqs, memProds);
    unsigned long long *d_pauliTargBitseqs;     cudaMalloc(&d_pauliTargBitseqs, memProds);
    unsigned long long *d_outcomeTargBitseqs;   cudaMalloc(&d_outcomeTargBitseqs, memProds);
    unsigned long long *d_baseBitseqs;          cudaMalloc(&d_baseBitseqs, memSamps);
    unsigned long long *d_outcomeBitseqs;       cudaMalloc(&d_outcomeBitseqs, memSamps);
    
    int bs = 128; // blocksize for GPU parallelisation
    
    // perform validation 
    int invalid = 0;
    int *d_invalid;
    cudaMalloc(&d_invalid, sizeof(int));
    cudaMemcpy(d_invalid, &invalid, sizeof(int), cudaMemcpyHostToDevice); 
    local_validateShadowPaulisKernel<<<ceil(numTotalPaulis/(qreal)bs), bs>>>(
        d_invalid, numQb, numTotalPaulis, d_pauliCodes, d_pauliTargs);
    local_validateShadowSampsKernel<<<ceil(numSamples*numQb/(qreal)bs), bs>>>(
        d_invalid, numSamples*numQb, d_sampleBases, d_sampleOutcomes);
    cudaMemcpy(&invalid, d_invalid, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_invalid);
    
    if (!invalid) {
        
        // prepare bit sequences
        local_prepareShadowPauliBitseqsKernel<<<ceil(numProds/(qreal)bs), bs>>>(
            d_pauliBitseqs, d_pauliTargBitseqs, d_outcomeTargBitseqs,
            numProds, d_numPaulisPerProd, d_pauliIndOffset, d_pauliCodes, d_pauliTargs);
        local_prepareShadowSampleBitseqsKernel<<<ceil(numSamples/(qreal)bs), bs>>>(
            d_baseBitseqs, d_outcomeBitseqs,
            numQb, numSamples, d_sampleBases, d_sampleOutcomes);
            
        // evaluate shadow expec values
        extension_calcExpecPauliProdsFromClassicalShadowKernel<<<ceil(numSamples/(qreal)bs), bs>>>(
            d_prodExpecVals,
            numProds, numSamples, numBatches, 
            d_numPaulisPerProd, d_baseBitseqs, d_pauliTargBitseqs, 
            d_pauliBitseqs, d_outcomeBitseqs, d_outcomeTargBitseqs);
            
        // copy expec vals back to RAM 
        cudaMemcpy(prodExpecVals.data(), d_prodExpecVals, memExpecVals, cudaMemcpyDeviceToHost);
    }

    // clean-up
    cudaFree(d_prodExpecVals);
    cudaFree(d_sampleOutcomes);
    cudaFree(d_pauliCodes);
    cudaFree(d_pauliTargs);
    cudaFree(d_numPaulisPerProd);
    cudaFree(d_pauliIndOffset);
    cudaFree(d_pauliBitseqs);
    cudaFree(d_pauliTargBitseqs);
    cudaFree(d_outcomeTargBitseqs);
    cudaFree(d_baseBitseqs);
    cudaFree(d_outcomeBitseqs);
    
    if (invalid)
        throw QuESTException("", "The input classical shadow, or the Pauli products, were invalid.");
}

