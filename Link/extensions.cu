/** @file
 * Contains GPU specific high-performance backend facilities,
 * which can be considered as extensions to the QuEST API.
 *
 * @author Tyson Jones
 */

#include "QuEST.h"
#include "QuEST_validation.h"



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
