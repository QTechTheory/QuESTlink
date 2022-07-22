/** @file 
 * Contains functions for processing analytic derivatives of circuits and channels.
 *
 * @author Tyson Jones
 */

#include "wstp.h"
#include "QuEST.h"
#include "QuEST_complex.h"
#include "QuEST_internal.h"
#include "QuEST_validation.h"

#include "errors.hpp"
#include "decoders.hpp"
#include "extensions.hpp"
#include "circuits.hpp"
#include "derivatives.hpp"
#include "link.hpp"



/*
 * Kraus map derivative superoperator preparation
 */

#define macro_populateKrausOperatorDeriv(superOp, ops, opDerivs, numOps, opDim) \
    /* clear the superop */ \
    for (int r=0; r < (opDim)*(opDim); r++) \
        for (int c=0; c < (opDim)*(opDim); c++) { \
            superOp->real[r][c] = 0; \
            superOp->imag[r][c] = 0; \
        } \
    /* add each op's contribution to the superop */ \
    for (int n=0; n<(numOps); n++) \
        /* superop += conjugate(op) (x) deriv(op), where (x) is a tensor product */ \
        for (int i = 0; i < (opDim); i++) \
            for (int j = 0; j < (opDim); j++) \
    			for (int k = 0; k < (opDim); k++) \
                    for (int l = 0; l < (opDim); l++) { \
                        superOp->real[i*(opDim) + k][j*(opDim) + l] += \
                              ops[n].real[i][j] * opDerivs[n].real[k][l] \
                            + ops[n].imag[i][j] * opDerivs[n].imag[k][l]; \
    					superOp->imag[i*(opDim) + k][j*(opDim) + l] += \
                              ops[n].real[i][j]*opDerivs[n].imag[k][l] \
                            - ops[n].imag[i][j]*opDerivs[n].real[k][l]; \
                    }            

void populateKrausSuperOperatorDerivN(ComplexMatrixN* superOp, ComplexMatrixN* ops, ComplexMatrixN* opDerivs, int numOps) {
    int opDim = 1 << ops[0].numQubits;
    macro_populateKrausOperatorDeriv(superOp, ops, opDerivs, numOps, opDim);
}



/*
 * decoherence derivativess (density matrix only)
 */

void local_mixTwoQubitDephasingDeriv(Qureg qureg, int targ1, int targ2, qreal probDeriv, Qureg workspace) {
    
    // TODO: remove workspace via bespoke extensions implementation

    cloneQureg(workspace, qureg);

    mixTwoQubitDephasing(qureg, targ1, targ2, probDeriv/2); // throws

    // qureg -> qureg - workspace
    Complex zero = {.real=0, .imag=0};
    Complex one = {.real=1, .imag=0};
    Complex neg = {.real=-1, .imag=0};
    setWeightedQureg(neg, workspace, zero, workspace, one, qureg);

    extension_addAdjointToSelf(qureg);
}

void local_mixDepolarisingDeriv(Qureg qureg, int targ, qreal probDeriv, Qureg workspace) {

    // TODO: remove workspace via bespoke extensions implementation

    cloneQureg(workspace, qureg);

    mixDepolarising(qureg, targ, probDeriv/2); // throws

    // qureg -> qureg - workspace
    Complex zero = {.real=0, .imag=0};
    Complex one = {.real=1, .imag=0};
    Complex neg = {.real=-1, .imag=0};
    setWeightedQureg(neg, workspace, zero, workspace, one, qureg);

    extension_addAdjointToSelf(qureg);
}

void local_mixTwoQubitDepolarisingDeriv(Qureg qureg, int targ1, int targ2, qreal probDeriv, Qureg workspace) {

    // TODO: remove workspace via bespoke extensions implementation

    cloneQureg(workspace, qureg);

    mixTwoQubitDepolarising(qureg, targ1, targ2, probDeriv/2); // throws

    // qureg -> qureg - workspace
    Complex zero = {.real=0, .imag=0};
    Complex one = {.real=1, .imag=0};
    Complex neg = {.real=-1, .imag=0};
    setWeightedQureg(neg, workspace, zero, workspace, one, qureg);

    extension_addAdjointToSelf(qureg);
}

void local_mixMultiQubitKrausMapDeriv(
    Qureg qureg, int* targs, int numTargs, 
    ComplexMatrixN* ops, ComplexMatrixN* opDerivs, int numOps
) {
    // validation must be done explicitly, since backend is called directly to avoid dens-matr full ops
    validateMultiTargets(qureg, targs, numTargs, "applyMultiQubitKrausSuperoperator");
    validateMultiQubitMatrixFitsInNode(qureg, 2*numTargs, "applyMultiQubitKrausSuperoperator");
    for (int i=0; i<numOps; i++) {
        validateMatrixInit(ops[i], "applyMultiQubitKrausSuperoperator");
        validateMatrixInit(opDerivs[i], "applyMultiQubitKrausSuperoperator");
        if (ops[i].numQubits != numTargs)
            throw QuESTException("", 
                "All Kraus matrices must have the same dimension, equal to 2^numTargs."); // throws
        // derivOps[n] gauaranteed equal sized to ops
    }
    
    ComplexMatrixN superOp = createComplexMatrixN(2*numTargs);
    populateKrausSuperOperatorDerivN(&superOp, ops, opDerivs, numOps);
    
    long long int ctrlMask = 0;
    int allTargs[MAX_NUM_TARGS_CTRLS]; // 2*numTargs
    for (int t=0; t < numTargs; t++) {
        allTargs[t] = targs[t];
        allTargs[t+numTargs] = targs[t] + qureg.numQubitsRepresented;
    }
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, allTargs, 2*numTargs, superOp);    
    
    destroyComplexMatrixN(superOp);
}



/*
 * gate derivatives (statevector & density matrix agnostic)
 */

void local_multiControlledMultiRotatePauliDeriv(
    Qureg qureg, 
    int* ctrls, int numCtrls, 
    int* targs, pauliOpType* paulis, int numTargs,
    qreal arg, qreal argDeriv
) {

    /* d/dx R[f(x), paulis] = (-i/2 * df/dx) [paulis] R[f(x), paulis] 
     * C_c[...] =  (...)|c=1><c=1|
     * note the R is a full operation, but the [paulis] and the projector act only 
     * on the lower qubits (of a density matrix).
     */

    // this call performs all input validation
    if (numCtrls > 0)
        multiControlledMultiRotatePauli(qureg, ctrls, numCtrls, targs, paulis, numTargs, arg); // throws
    else
        multiRotatePauli(qureg, targs, paulis, numTargs, arg); // throws

    for (int c=0; c<numCtrls; c++)
        statevec_collapseToKnownProbOutcome(qureg, ctrls[c], 1, 1);

    for (int t=0; t < numTargs; t++) {
        if (paulis[t] == PAULI_X)
            statevec_pauliX(qureg, targs[t]);
        if (paulis[t] == PAULI_Y)
            statevec_pauliY(qureg, targs[t]);
        if (paulis[t] == PAULI_Z)
            statevec_pauliZ(qureg, targs[t]);
    }

    extension_applyImagFactor(qureg, -argDeriv/2);

    if (qureg.isDensityMatrix)
        extension_addAdjointToSelf(qureg);
}

void local_multiControlledPhaseShiftDeriv(
    Qureg qureg, int* ctrls, int numCtrls, qreal arg, qreal argDeriv
) {

    /* d/dx C_c[Ph_t[f(x)] = d/dx Ph_{c & t}[f(x)] = (i df/dx) |c,t=1><c,t=1| Ph_{c & t}[f(x)]
     * where Ph is a full operation, but the projector applies only to the lower qubits 
     * of a density matrix.
     */

    // this call performs all input validation
    multiControlledPhaseShift(qureg, ctrls, numCtrls, arg); // throws

    for (int c=0; c<numCtrls; c++)
        statevec_collapseToKnownProbOutcome(qureg, ctrls[c], 1, 1);

    extension_applyImagFactor(qureg, argDeriv);

    if (qureg.isDensityMatrix)
        extension_addAdjointToSelf(qureg);
}

void local_multiControlledMultiQubitMatrixDeriv(
    Qureg qureg, 
    int* ctrls, int numCtrls, 
    int* targs, int numTargs,
    ComplexMatrixN matr, ComplexMatrixN matrDeriv
) {
    // validation must be done explicitly, since backend is called directly to avoid dens-matr full ops
    validateMultiQubitMatrix(qureg, matrDeriv, numTargs, "applyGateMatrixN");
    if (qureg.isDensityMatrix)
        validateMultiQubitMatrix(qureg, matr, numTargs, "applyGateMatrixN"); // only featured in dens matr picture
    if (numCtrls > 0)
        validateMultiControlsMultiTargets(qureg, ctrls, numCtrls, targs, numTargs, "applyMultiControlledGateMatrixN");
    else
        validateMultiTargets(qureg, targs, numTargs, "applyGateMatrixN");

    // state-vector (low region of density matrix)
    long long int ctrlMask = getQubitBitMask(ctrls, numCtrls);
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targs, numTargs, matrDeriv);

    for (int c=0; c<numCtrls; c++)
        statevec_collapseToKnownProbOutcome(qureg, ctrls[c], 1, 1);

    // density matrix (upper region)
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(matr);
        statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask<<shift, targs, numTargs, matr);
        setConjugateMatrixN(matr);
        shiftIndices(targs, numTargs, - shift);
    }

    if (qureg.isDensityMatrix)
        extension_addAdjointToSelf(qureg);
}

void local_globalPhaseDeriv(Qureg qureg, qreal phaseDeriv) {
    
    if (qureg.isDensityMatrix)
        initBlankState(qureg);
    else
        extension_applyImagFactor(qureg, phaseDeriv);
}



/*
 * Gate methods
 */
 
pauliOpType* local_preparePauliCache(pauliOpType pauli, int numPaulis) {
    
    static pauliOpType paulis[MAX_NUM_TARGS_CTRLS]; 
    for (int i=0; i < numPaulis; i++)
        paulis[i] = pauli;
    return paulis;
}

void Gate::applyDerivTo(Qureg qureg, qreal* derivParams, int numDerivParams) {
    
    // TODO: remove workspace via bespoke extensions implementation
    // (especially because it gets leaked upon error throw)
    Qureg workspace = createCloneQureg(qureg, env);
    
    switch(opcode) {
            
        case OPCODE_Rx : {
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Rx", numDerivParams, 1); // throws
            local_multiControlledMultiRotatePauliDeriv(
                qureg, ctrls, numCtrls, 
                targs, local_preparePauliCache(PAULI_X, numTargs), numTargs,
                params[0], derivParams[0]);
        }
            break;
            
        case OPCODE_Ry : {
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Ry", numDerivParams, 1); // throws
            local_multiControlledMultiRotatePauliDeriv(
                qureg, ctrls, numCtrls, 
                targs, local_preparePauliCache(PAULI_Y, numTargs), numTargs,
                params[0], derivParams[0]);
        }
            break;
            
        case OPCODE_Rz : {
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Rz", numDerivParams, 1); // throws
            local_multiControlledMultiRotatePauliDeriv(
                qureg, ctrls, numCtrls, 
                targs, local_preparePauliCache(PAULI_Z, numTargs), numTargs,
                params[0], derivParams[0]);
        }
            break;
            
        case OPCODE_R: {
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("R", numDerivParams, 1); // throws
                
            pauliOpType paulis[MAX_NUM_TARGS_CTRLS]; // numTargs
            for (int p=0; p < numTargs; p++)
                paulis[p] = (pauliOpType) ((int) params[1+p]);
                
            local_multiControlledMultiRotatePauliDeriv(
                qureg, ctrls, numCtrls, 
                targs, paulis, numTargs,
                params[0], derivParams[0]);
        }
            break;
        
        case OPCODE_U :     // intentional fallthrough
        case OPCODE_UNonNorm : 
        case OPCODE_Matr : { ; 
            int reqNumDerivParams = local_getNumScalarsToFormMatrix(numTargs);
            if (numDerivParams != reqNumDerivParams)
                throw local_wrongNumDerivParamsExcep("U", numDerivParams, reqNumDerivParams); // throws
            
            ComplexMatrixN matr = createComplexMatrixN(numTargs);
            ComplexMatrixN matrDeriv = createComplexMatrixN(numTargs);
            local_setMatrixNFromFlatList(params, matr, numTargs);
            local_setMatrixNFromFlatList(derivParams, matrDeriv, numTargs);
            
            local_multiControlledMultiQubitMatrixDeriv(
                qureg, ctrls, numCtrls, targs, numTargs, matr, matrDeriv);
            
            destroyComplexMatrixN(matr);
            destroyComplexMatrixN(matrDeriv);
        }
            break;
            
        case OPCODE_Deph :
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Dephasing", numDerivParams, 1); // throws
            if (numTargs == 1)
                extension_mixDephasingDeriv(qureg, targs[0], derivParams[0]);
            if (numTargs == 2)
                local_mixTwoQubitDephasingDeriv(qureg, targs[0], targs[1], derivParams[0], workspace);
            break;
            
        case OPCODE_Depol :
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Depolarising", numDerivParams, 1); // throws
            if (numTargs == 1)
                local_mixDepolarisingDeriv(qureg, targs[0], derivParams[0], workspace);
            if (numTargs == 2)
                local_mixTwoQubitDepolarisingDeriv(qureg, targs[0], targs[1], derivParams[0], workspace);
            break;
            
        case OPCODE_Damp :
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Damping", numDerivParams, 1); // throws
            extension_mixDampingDeriv(qureg, targs[0], params[0], derivParams[0]);
            break;
            
        case OPCODE_Kraus :         // intentional fallthrough
        case OPCODE_KrausNonTP : { ;
            int numOps = (int) params[0];
            int reqNumDerivParams = numOps * local_getNumScalarsToFormMatrix(numTargs);
            if (numDerivParams != reqNumDerivParams)
                throw local_wrongNumDerivParamsExcep("Kraus map", numDerivParams, reqNumDerivParams); // throws
            
            ComplexMatrixN* ops = (ComplexMatrixN*) malloc(numOps * sizeof *ops);
            ComplexMatrixN* opDerivs = (ComplexMatrixN*) malloc(numOps * sizeof *opDerivs);
            local_createManyMatrixNFromFlatList(&params[1], ops, numOps, numTargs);
            local_createManyMatrixNFromFlatList(&derivParams[1], opDerivs, numOps, numTargs);
            
            local_mixMultiQubitKrausMapDeriv(qureg, targs, numTargs, ops, opDerivs, numOps);
            
            for (int i=0; i<numOps; i++) {
                destroyComplexMatrixN(ops[i]);
                destroyComplexMatrixN(opDerivs[i]);
            }
            free(ops);
            free(opDerivs);
        }
            break;
            
        case OPCODE_G :
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Global phase", numDerivParams, 1); // throws
            local_globalPhaseDeriv(qureg, derivParams[0]);
            break;
            
        case OPCODE_Ph : {
            if (numDerivParams != 1)
                throw local_wrongNumDerivParamsExcep("Phase gate", numDerivParams, 1); // throws
                
            int* qubitCache = local_prepareCtrlCache(ctrls, numCtrls, -1);
            for (int i=0; i<numTargs; i++)
                qubitCache[numCtrls+i] = targs[i];
            
            local_multiControlledPhaseShiftDeriv(
                qureg, qubitCache, numCtrls+numTargs, params[0], derivParams[0]);
        }
            break;
            
        default:            
            throw QuESTException("", "The circuit contained a parameterised gate which does not "
                "have a known analytic derivative."); // throws
    }

    destroyQureg(workspace, env);
}



/*
 * DerivTerm methods 
 */

void DerivTerm::init(Gate gate, int gateInd, int varInd, qreal* derivParams, int numDerivParams) {
    
    this->gate = gate;
    this->gateInd = gateInd;
    this->varInd = varInd;
    this->derivParams = derivParams;
    this->numDerivParams = numDerivParams;
}

void DerivTerm::applyTo(Qureg qureg) {
    
    gate.applyDerivTo(qureg, derivParams, numDerivParams);
}



/*
 * DerivCircuit methods 
 */

void DerivCircuit::applyTo(Qureg* quregs, int numQuregs, Qureg initQureg) {
    
    /* A reader may desire to optimise this function by cloning partial states 
     * between the derivQuregs so that no single gate is performed more 
     * than once per local-diff-variable. This would only offer a factor 2x 
     * speedup (since the remaining gates of each partial state must be sequentially
     * applied), but is anyway precluded by our support of repeated parameters / 
     * chain rule through use of a workspace. Pity!
     */
    
    for (int q=0; q<numQuregs; q++)
        initBlankState(quregs[q]);
    
    Qureg workspace = createCloneQureg(initQureg, env);
    
    for (int t=0; t<numTerms; t++) {
        
        DerivTerm term = terms[t];
        int gateInd = term.getGateInd();
        
        cloneQureg(workspace, initQureg);
        circuit->applySubTo(workspace, 0, gateInd); // throws
        term.applyTo(workspace); // throws
        circuit->applySubTo(workspace, gateInd+1, circuit->getNumGates()); // throws
            
        Qureg qureg = quregs[term.getVarInd()];
        Complex zero = {.real=0, .imag=0};
        Complex one = {.real=1, .imag=0};    
        setWeightedQureg(zero, workspace, one, workspace, one, qureg);
    }

    // clean up
    destroyQureg(workspace, env);
}

void DerivCircuit::calcDerivEnergiesStateVec(qreal* eneryGrad, PauliHamil hamil, Qureg initQureg) {
    
    if (!circuit->isUnitary())
        throw QuESTException("", "The given circuit must be composed strictly of unitary gates "
            "(and ergo exclude gates like Matr[] and P[]), in order to return a valid real "
            "observable gradient. Please instead use CalcQuregDerivs[]"); // throws
    
    // prepare lambda = H circuit(init) and phi=circuit(init)
    Qureg workLambda = createCloneQureg(initQureg, env);
    circuit->applyTo(workLambda);
    Qureg workPhi = createCloneQureg(workLambda, env);
    applyPauliHamil(workPhi, hamil, workLambda);
    
    // prepare workMu in arbitrary state
    Qureg workMu = createQureg(initQureg.numQubitsRepresented, env);

    // clear energies
    for (int i=0; i<numVars; i++)
        eneryGrad[i] = 0;
        
    // explicitly catch and rethrow errors to force clean-up of above quregs
    try {
        
        for (int t=numTerms-1; t>=0; t--) {
            
            DerivTerm derivTerm = terms[t];
            int gateInd = derivTerm.getGateInd();
            int varInd = derivTerm.getVarInd();
            
            // remove all gates >= gateInd (not removed by previous iteration) from workPhi
            circuit->applyDaggerSubTo(workPhi, gateInd, 
                (t<numTerms-1)? terms[t+1].getGateInd() : circuit->getNumGates());

            // add all (daggered) gates > gateInd (not added by previous iteration) to workLambda
            circuit->applyDaggerSubTo(workLambda, gateInd + 1,
                (t<numTerms-1)? terms[t+1].getGateInd() + 1 : circuit->getNumGates());
                
            cloneQureg(workMu, workPhi);
            derivTerm.applyTo(workMu);
            eneryGrad[varInd] += 2 * calcInnerProduct(workLambda, workMu).real;        
        }
        
    } catch (QuESTException& err) {
        
        // clean-up and rethrow
        destroyQureg(workLambda, env);
        destroyQureg(workPhi, env);
        destroyQureg(workMu, env);
        throw;
    }
    
    // clean-up
    destroyQureg(workLambda, env);
    destroyQureg(workPhi, env);
    destroyQureg(workMu, env);
}

void DerivCircuit::calcDerivEnergiesDensMatr(qreal* energyGrad, PauliHamil hamil, Qureg initQureg) {
    
    // initQureg may be a statevector or a density matrix
    Qureg qureg = createDensityQureg(initQureg.numQubitsRepresented, env);
    Qureg workspace = createCloneQureg(qureg, env);
    
    // clear energies
    for (int i=0; i<numVars; i++)
        energyGrad[i] = 0;
    
    // explicitly catch and rethrow errors to force clean-up of above quregs
    try {
        
        // iterate each differential term (those produced after chain-rule expansion)
        for (int t=0; t<numTerms; t++) {
            
            DerivTerm term = terms[t];
            int gateInd = term.getGateInd();
            int varInd = term.getVarInd();
            
            // set qureg = initQureg
            if (initQureg.isDensityMatrix)
                cloneQureg(qureg, initQureg);
            else
                initPureState(qureg, initQureg);
            
            // produce a single term of (d circuit(qureg) / d var)
            circuit->applySubTo(qureg, 0, gateInd); // throws
            term.applyTo(qureg); // throws
            circuit->applySubTo(qureg, gateInd+1, circuit->getNumGates()); // throws
            
            // add this term to (d circuit(qureg) / d var)
            energyGrad[varInd] += calcExpecPauliHamil(qureg, hamil, workspace);
        }
        
    } catch (QuESTException& err) {
        
        // clean-up and rethrow
        destroyQureg(qureg, env);
        destroyQureg(workspace, env);
        throw;
    }

    // clean-up
    destroyQureg(qureg, env);
    destroyQureg(workspace, env);
}

void DerivCircuit::calcDerivEnergies(qreal* energyGrad, PauliHamil hamil, Qureg initQureg, bool isPureCirc) {
    
    if (isPureCirc && !initQureg.isDensityMatrix)
        calcDerivEnergiesStateVec(energyGrad, hamil, initQureg);
    else
        calcDerivEnergiesDensMatr(energyGrad, hamil, initQureg);
}

void DerivCircuit::calcGeometricTensorStateVec(qcomp** tensor, Qureg initQureg) {
    
    Qureg quregDiag = createCloneQureg(initQureg, env);
    Qureg quregSuffix = createCloneQureg(initQureg, env);
    Qureg quregPrefix = createCloneQureg(initQureg, env);
    Qureg quregDeriv = createCloneQureg(initQureg, env);
    
    // clear tensor
    for (int i=0; i<numVars; i++)
        for (int j=0; j<numVars; j++)
            tensor[i][j] = 0;
            
    // separate array for the Berry connection terms
    qcomp* berries = (qcomp*) calloc(numVars, sizeof* berries);
    
    try {
        int indOfLastGateOnDiag = -1;
        
        for (int t=0; t<numTerms; t++) {
            
            DerivTerm rowDerivTerm = terms[t];
            int rowGateInd = rowDerivTerm.getGateInd();
            int rowVarInd = rowDerivTerm.getVarInd();

            // apply all prior gates to |diag>, such that 
            // |diag> = U_(r-1) ... U2 U1 |in>
            circuit->applySubTo(quregDiag, indOfLastGateOnDiag+1, rowGateInd); // throws
            indOfLastGateOnDiag = rowGateInd-1;
        
            // |deriv> = (dU_r/dx) U_(r-1) ... U2 U1 |in>
            cloneQureg(quregDeriv, quregDiag);
            rowDerivTerm.applyTo(quregDeriv); // throws
            
            // <deriv|deriv> = || (dU_r/dx) U_(r-1) ... U2 U1 |in> ||^2
            Complex norm = calcInnerProduct(quregDeriv, quregDeriv);
            tensor[rowVarInd][rowVarInd] += fromComplex(norm);
            
            // <prefix| = <in| U1^ U2^ ... U_(r-1)^ (dU_r/dx)^
            cloneQureg(quregPrefix, quregDeriv);
            
            // |suffix> = U_(r-1) ... U2 U1 |in>
            cloneQureg(quregSuffix, quregDiag);
            
            // aside: compute single Berry connection term (safely modify quregDeriv)
            cloneQureg(quregDeriv, quregDiag);
            circuit->applySubTo(quregDeriv, indOfLastGateOnDiag+1, rowGateInd+1);
            Complex berry = calcInnerProduct(quregPrefix, quregDeriv);
            berries[rowVarInd] += fromComplex(berry);
            
            // mark that prefix needs only gates before U_r (inclusive)
            int indOfLastGateOnPrefix = rowGateInd+1;   // as added
            int indOfLastGateOnSuffix = rowGateInd-1;   // as remaining
            
            for (int s=t-1; s>=0; s--) {
                
                DerivTerm colDerivTerm = terms[s];
                int colGateInd = colDerivTerm.getGateInd();
                int colVarInd = colDerivTerm.getVarInd();

                // <prefix| = <in| U1^ U2^ ... U_(r-1)^ (dU_r/dx)^ U_r U_(r-1) ... U_(c+1)                    
                circuit->applyDaggerSubTo(quregPrefix, colGateInd+1, indOfLastGateOnPrefix);
                indOfLastGateOnPrefix = colGateInd + 1;
                
                // |suffix> = U_(c-1) ... U2 U1 |in>
                circuit->applyInverseSubTo(quregSuffix, colGateInd, indOfLastGateOnSuffix+1);
                indOfLastGateOnSuffix = colGateInd - 1;
                
                // |deriv> = (dU_c/dx) U_(c-1) ... U2 U1 |in>
                cloneQureg(quregDeriv, quregSuffix);
                colDerivTerm.applyTo(quregDeriv);

                // <prefix|deriv> = 
                // <in| U1^ U2^ ... U_(r-1)^ (dU_r/dx)^ U_r U_(r-1) ... U_(c+1) (dU_c/dx) U_(c-1) ... U2 U1 |in>
                Complex prod = calcInnerProduct(quregPrefix, quregDeriv);
                tensor[rowVarInd][colVarInd] += fromComplex(prod);
                tensor[colVarInd][rowVarInd] += conj(fromComplex(prod));
            }
        }
        
        // add Berry connections to the geometric tensor
        for (int r=0; r<numVars; r++)
            for (int c=0; c<numVars; c++)
                tensor[r][c] -= berries[r] * conj(berries[c]);
        
    } catch (QuESTException& err) {
        
        // clean-up and rethrow
        free(berries);
        destroyQureg(quregDiag, env);
        destroyQureg(quregSuffix, env);
        destroyQureg(quregPrefix, env);
        destroyQureg(quregDeriv, env);
        throw;
    }

    // clean-up
    free(berries);
    destroyQureg(quregDiag, env);
    destroyQureg(quregSuffix, env);
    destroyQureg(quregPrefix, env);
    destroyQureg(quregDeriv, env);
}

void DerivCircuit::calcGeometricTensorDensMatr(qcomp** tensor, Qureg initQureg) {
    
    /* TODO:
     * Implement quantum Fisher information matrix?
     * - Thm 1 https://arxiv.org/pdf/1801.00945.pdf gives QFI in terms of inv(rho) and d rho/dx
     * - inv(rho) could be estimated with https://arxiv.org/pdf/cs/0412107.pdf
     */
    
    throw QuESTException("", "This facility is not yet available for channels or density matrices.");
}

void DerivCircuit::calcGeometricTensor(qcomp** tensor, Qureg initQureg, bool isPureCirc) {
    
    if (isPureCirc &&  !initQureg.isDensityMatrix)
        calcGeometricTensorStateVec(tensor, initQureg); // throws
    else
        calcGeometricTensorDensMatr(tensor, initQureg); // throws
}

DerivCircuit::~DerivCircuit() {
    
    freeMMA();
    delete circuit;
    delete[] terms;
}



/*
 * interfacing
 */

void internal_calcQuregDerivs(int initQuregId) {
    
    // get qureg ids (one for each var)
    int* quregIds;
    int numQuregs;
    WSGetInteger32List(stdlink, &quregIds, &numQuregs); // must free
    
    // load the circuit and derivative descriptions (local so no need to explicitly delete)
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA();
    
    // validate quregs (must do so after loading derivCircuit from MMA so those packets are flushed)
    try {
        local_throwExcepIfQuregNotCreated(initQuregId); // throws
        
        int numQb = quregs[initQuregId].numQubitsRepresented;
        int isDens = quregs[initQuregId].isDensityMatrix;
        for (int q=0; q<numQuregs; q++) {
            local_throwExcepIfQuregNotCreated(quregIds[q]); // throws
            if (quregIds[q] == initQuregId)
                throw QuESTException("", "The initial state qureg must not also be one of the derivative quregs"); // throws 
            if (quregs[quregIds[q]].numQubitsRepresented != numQb)
                throw QuESTException("", "All derivative quregs must have the same dimension as the initial state."); // throws
            if (quregs[quregIds[q]].isDensityMatrix != isDens)
                throw QuESTException("", "Quregs must be all state-vectors or all density-matrices"); // throws
        }
    } catch (QuESTException& err) {
        local_sendErrorAndFail("CalcQuregDerivs", err.message);
        return;
    }
    
    Qureg initQureg = quregs[initQuregId];
    
    Qureg* derivQuregs = (Qureg*) malloc(numQuregs * sizeof *derivQuregs);
    for (int q=0; q<numQuregs; q++)
        derivQuregs[q] = quregs[quregIds[q]];
        
    try {
        derivCirc.applyTo(derivQuregs, numQuregs, initQureg); // throws
        
        WSPutInteger(stdlink, initQuregId);
        
    } catch (QuESTException& err) {
        
        // report error, depending on type
        if (err.thrower == "")
            local_sendErrorAndFail("CalcQuregDerivs", err.message);
        else if (err.thrower == "Abort")
            local_sendErrorAndAbort("CalcQuregDerivs", err.message);
        else 
            local_sendErrorAndFail("CalcQuregDerivs", 
                "Problem in " + err.thrower + ": " + err.message);
    } 
    
    // clean-up even in event of error 
    free(derivQuregs);
}

void internal_calcExpecPauliSumDerivs(int initQuregId, int isPureCirc) {
    
    // load the circuit and deriv spec from MMA
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA(); // local, so desconstructor automatic
    
    // load Hamiltonian from MMA (and also validate initQuregId)
    PauliHamil hamil;
    try {
        hamil = local_loadPauliHamilForQuregFromMMA(initQuregId); // throws
        
    } catch (QuESTException& err) {
        
        local_sendErrorAndFail("CalcExpecPauliSumDerivs", err.message);
        return;
    }
    
    // prepare energy grad vector (malloc onto heap to avoid stack size limits)
    int numDerivs = derivCirc.getNumVars();
    qreal* energyGrad = (qreal*) malloc(numDerivs * sizeof *energyGrad);
    
    try {    
        derivCirc.calcDerivEnergies(energyGrad, hamil, quregs[initQuregId], isPureCirc); // throws
        
        WSPutReal64List(stdlink, energyGrad, numDerivs);
        
    } catch (QuESTException& err) {
        
        local_sendErrorAndFail("CalcExpecPauliSumDerivs", err.message);
    }

    // clean-up even despite errors
    free(energyGrad);
    local_freePauliHamil(hamil);
}

void internal_calcGeometricTensor(int initQuregId, int isPureCirc) {
    
    // load the circuit and deriv spec from MMA
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA(); // local, so desconstructor automatic
    
    int tensorDim = derivCirc.getNumVars();
    qcomp** tensor = (qcomp**) malloc(tensorDim * sizeof *tensor);
    for (int i=0; i<tensorDim; i++)
        tensor[i] = (qcomp*) malloc(tensorDim * sizeof **tensor);
    
    try {
        
        local_throwExcepIfQuregNotCreated(initQuregId); // throws

        derivCirc.calcGeometricTensor(tensor, quregs[initQuregId], isPureCirc); // throws
        
        local_sendMatrixToMMA(tensor, tensorDim);
    
    } catch (QuESTException& err) {
        
        local_sendErrorAndFail("CalcGeometricTensor", err.message);
    }
    
    // clean-up, even if error
    for (int i=0; i<tensorDim; i++)
        free(tensor[i]);
    free(tensor);
}
