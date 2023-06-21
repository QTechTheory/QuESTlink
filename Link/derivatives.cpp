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
#include "utilities.hpp"
#include "derivatives.hpp"
#include "link.hpp"



/*
 * Kraus derivatives (in addition to those in extensions.cpp)
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

void populateKrausSuperOperatorDerivN(
    ComplexMatrixN* superOp, ComplexMatrixN* ops, ComplexMatrixN* opDerivs, int numOps
) {
    int opDim = 1 << ops[0].numQubits;
    macro_populateKrausOperatorDeriv(superOp, ops, opDerivs, numOps, opDim);
}

void local_mixMultiQubitKrausMapDeriv(
    Qureg qureg, int* targs, int numTargs, 
    ComplexMatrixN* ops, ComplexMatrixN* opDerivs, int numOps
) {
    // validation must be done explicitly, since backend is called directly to avoid dens-matr full ops
    validateMultiTargets(qureg, targs, numTargs, "applyMultiQubitKrausSuperoperator"); // throws
    validateMultiQubitMatrixFitsInNode(qureg, 2*numTargs, "applyMultiQubitKrausSuperoperator"); // throws
    for (int i=0; i<numOps; i++) {
        validateMatrixInit(ops[i], "applyMultiQubitKrausSuperoperator"); // throws
        validateMatrixInit(opDerivs[i], "applyMultiQubitKrausSuperoperator"); // throws
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
    // d/dx R[f(x), paulis] = (-i/2 * df/dx) [paulis] R[f(x), paulis] 
    // C_c[...] =  (...)|c=1><c=1|
    // note the R is a full operation, but the [paulis] and the projector act only 
    // on the lower qubits (of a density matrix).

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
    // d/dx C_c[Ph_t[f(x)] = d/dx Ph_{c & t}[f(x)] = (i df/dx) |c,t=1><c,t=1| Ph_{c & t}[f(x)]
    // where Ph is a full operation, but the projector applies only to the lower qubits 
    // of a density matrix.

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
    ComplexMatrixN matr, ComplexMatrixN matrDeriv,
    bool leftApplyOnly
) {    
    // all validation of ctrls, targs and matrDeriv is done by apply*MatrixN call below, but
    // validation of matr (only used by density mtrix) must be done explicitly since we call backend func directly
    if (qureg.isDensityMatrix)
        validateMultiQubitMatrix(qureg, matr, numTargs, "applyGateMatrixN"); // throws

    // |psi> -> C[D(M)]|psi>,  |rho> -> C[D(M)]|rho>   (notice controls C is outside of deriv D)
    if (numCtrls > 0)
        applyMultiControlledMatrixN(qureg, ctrls, numCtrls, targs, numTargs, matrDeriv); // throws
    else 
        applyMatrixN(qureg, targs, numTargs, matrDeriv); // throws

    // |psi> -> D(C[M])|psi>,  |rho> -> D(C[M])|rho>
    for (int c=0; c<numCtrls; c++)
        statevec_collapseToKnownProbOutcome(qureg, ctrls[c], 1, 1);
        
    // subsequent Jamiolkowski-isomorphism and chain rule don't occur if only left-applying D[gate], like for Matr[]
    if (leftApplyOnly)
        return;

    // |rho> -> conj(C[M]) (x) D(C[M])|rho>
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        long long int ctrlMask = getQubitBitMask(ctrls, numCtrls);
        
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(matr);
        statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask<<shift, targs, numTargs, matr);
        setConjugateMatrixN(matr);
        shiftIndices(targs, numTargs, - shift);
    }

    // |rho> -> conj(C[M]) (x) D(C[M]) |rho> + C[M] (x) conj(D(C[M])) |rho>
    if (qureg.isDensityMatrix)
        extension_addAdjointToSelf(qureg);
}

void local_subDiagonalOpDeriv(
    Qureg qureg, int* targs, int numTargs, 
    SubDiagonalOp op, SubDiagonalOp opDeriv,
    bool leftApplyOnly
) {
    // |psi> -> D(op)|psi>,  |rho> -> D(op)|rho>
    applySubDiagonalOp(qureg, targs, numTargs, opDeriv);
    
    // subsequent Jamiolkowski-isomorphism and chain rule don't occur if only left-applying D[gate], like for Matr[]
    if (leftApplyOnly)
        return;
    
    // |rho> -> conj(op) (x) D(op) |rho>
    if (qureg.isDensityMatrix) {
        int conj = 1;
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, op.numQubits, shift);
        statevec_applySubDiagonalOp(qureg, targs, op, conj);
        shiftIndices(targs, op.numQubits, - shift);
    }
    
    // |rho> -> conj(op) (x) D(op) |rho> + op (x) conj(D(op)) |rho>
    if (qureg.isDensityMatrix)
        extension_addAdjointToSelf(qureg);
}

void local_globalPhaseDeriv(Qureg qureg, qreal phaseDeriv) {
    
    if (qureg.isDensityMatrix)
        initBlankState(qureg);
    else
        extension_applyImagFactor(qureg, phaseDeriv);
}

void local_factorDeriv(Qureg qureg, qreal facDerivRe, qreal facDerivIm) {
    
    Complex zero;  zero.real = 0;           zero.imag = 0;
    Complex fac;    fac.real = facDerivRe;   fac.imag = facDerivIm;
    setWeightedQureg(zero, qureg, zero, qureg, fac, qureg);
}



/*
 * Gate methods
 */

void Gate::applyDerivTo(Qureg qureg, qreal* derivParams, int numDerivParams) {
    
    validate(); // throws (and injects gate::getSyntax() into exception.thrower)
    
    // catch and rethrow errors so that we can inject the gate's Mathematica syntax 
    // into the exception, regardless of whether the exception is thrown by explicit 
    // validation or the QuEST backend. We ergo do not bother setting gateSyntax in the 
    // explicit validation throws below, replacing them with emoty strings.
    try {

        switch(opcode) {
                
            case OPCODE_Rx : {
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                local_multiControlledMultiRotatePauliDeriv(
                    qureg, ctrls, numCtrls, targs, 
                    local_preparePauliCache(PAULI_X, numTargs), numTargs,
                    params[0], derivParams[0]); // throws
            }
                break;
                
            case OPCODE_Ry : {
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                local_multiControlledMultiRotatePauliDeriv(
                    qureg, ctrls, numCtrls, targs, 
                    local_preparePauliCache(PAULI_Y, numTargs), numTargs,
                    params[0], derivParams[0]);  // throws
            }
                break;
                
            case OPCODE_Rz : {
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                local_multiControlledMultiRotatePauliDeriv(
                    qureg, ctrls, numCtrls, targs, 
                    local_preparePauliCache(PAULI_Z, numTargs), numTargs,
                    params[0], derivParams[0]);  // throws
            }
                break;
                
            case OPCODE_R:
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                local_multiControlledMultiRotatePauliDeriv(
                    qureg, ctrls, numCtrls, targs, 
                    local_preparePauliCache(&params[1], numTargs), numTargs,
                    params[0], derivParams[0]);  // throws
                break;
            
            case OPCODE_U :     // intentional fallthrough
            case OPCODE_UNonNorm : 
            case OPCODE_Matr : {
                bool leftApplyOnly = (opcode == OPCODE_Matr);
                if (local_isEncodedMatrix(params[0])) {
                    int reqNumDerivParams = local_getNumRealScalarsToFormMatrix(numTargs);
                    if (numDerivParams != reqNumDerivParams)
                        throw local_wrongNumDerivParamsExcep("", numDerivParams, reqNumDerivParams); // throws
                    
                    ComplexMatrixN matr = createComplexMatrixN(numTargs);
                    ComplexMatrixN matrDeriv = createComplexMatrixN(numTargs);
                    local_setMatrixNFromFlatList(&params[1], matr, numTargs);
                    local_setMatrixNFromFlatList(derivParams, matrDeriv, numTargs);
                    
                    local_multiControlledMultiQubitMatrixDeriv(
                        qureg, ctrls, numCtrls, targs, numTargs, matr, matrDeriv, leftApplyOnly); // throws
                    
                    destroyComplexMatrixN(matr);
                    destroyComplexMatrixN(matrDeriv);
                }
                if (local_isEncodedVector(params[0])) {
                    int reqNumDerivParams = local_getNumRealScalarsToFormDiagonalMatrix(numTargs);
                    if (numDerivParams != reqNumDerivParams)
                        throw local_wrongNumDerivParamsExcep("", numDerivParams, reqNumDerivParams); // throws
                        
                    SubDiagonalOp op = createSubDiagonalOp(numTargs);
                    SubDiagonalOp opDeriv = createSubDiagonalOp(numTargs);
                    local_setSubDiagonalOpFromFlatList(&params[1], op);
                    local_setSubDiagonalOpFromFlatList(derivParams, opDeriv);
                    
                    local_subDiagonalOpDeriv(qureg, targs, numTargs, op, opDeriv, leftApplyOnly);
                    
                    destroySubDiagonalOp(op);
                    destroySubDiagonalOp(opDeriv);
                }
            };
                break;
                
            case OPCODE_Deph :
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                if (numTargs == 1)
                    extension_mixDephasingDeriv(qureg, targs[0], derivParams[0]); // throws
                if (numTargs == 2)
                    extension_mixTwoQubitDephasingDeriv(qureg, targs[0], targs[1], derivParams[0]);  // throws
                break;
                
            case OPCODE_Depol :
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                if (numTargs == 1)
                    extension_mixDepolarisingDeriv(qureg, targs[0], derivParams[0]);  // throws
                if (numTargs == 2)
                    extension_mixTwoQubitDepolarisingDeriv(qureg, targs[0], targs[1], derivParams[0]);  // throws
                break;
                
            case OPCODE_Damp :
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                extension_mixDampingDeriv(qureg, targs[0], params[0], derivParams[0]);  // throws
                break;
                
            case OPCODE_Kraus :         // intentional fallthrough
            case OPCODE_KrausNonTP : { ;
                int numOps = (int) params[0];
                int reqNumDerivParams = numOps * local_getNumRealScalarsToFormMatrix(numTargs);
                if (numDerivParams != reqNumDerivParams)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, reqNumDerivParams); // throws
                
                ComplexMatrixN* ops = (ComplexMatrixN*) malloc(numOps * sizeof *ops);
                ComplexMatrixN* opDerivs = (ComplexMatrixN*) malloc(numOps * sizeof *opDerivs);
                local_createManyMatrixNFromFlatList(&params[1], ops, numOps, numTargs);
                local_createManyMatrixNFromFlatList(&derivParams[1], opDerivs, numOps, numTargs);
                
                local_mixMultiQubitKrausMapDeriv(qureg, targs, numTargs, ops, opDerivs, numOps);  // throws
                
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
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                local_globalPhaseDeriv(qureg, derivParams[0]);  // throws
                break;
                
            case OPCODE_Fac :
                if (numDerivParams != 2)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 2); // throws
                local_factorDeriv(qureg, derivParams[0], derivParams[1]);
                break;
                
            case OPCODE_Ph : {
                if (numDerivParams != 1)
                    throw local_wrongNumDerivParamsExcep("", numDerivParams, 1); // throws
                    
                int* qubitCache = local_prepareCtrlCache(ctrls, numCtrls, -1);
                for (int i=0; i<numTargs; i++)
                    qubitCache[numCtrls+i] = targs[i];
                
                local_multiControlledPhaseShiftDeriv(
                    qureg, qubitCache, numCtrls+numTargs, params[0], derivParams[0]);  // throws
            }
                break;
                
            default:            
                throw QuESTException("", 
                    "The operator family (" + getSymb() + ") has no known analytic derivative."); // throws, caught by below
        }
        
    } catch (QuESTException& err) {
        
        err.thrower = getSyntax();
        throw;
    }
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
    
    gate.applyDerivTo(qureg, derivParams, numDerivParams); // throws
}



/*
 * DerivCircuit methods 
 */

void DerivCircuit::applyTo(Qureg* quregs, int numQuregs, Qureg initQureg, Qureg workspace) {
    
    /* A reader may desire to optimise this function by cloning partial states 
     * between the derivQuregs so that no single gate is performed more 
     * than once per local-diff-variable. This would only offer a factor 2x 
     * speedup (since the remaining gates of each partial state must be sequentially
     * applied), but is anyway precluded by our support of repeated parameters / 
     * chain rule through use of a workspace. Pity!
     */
    
    for (int q=0; q<numQuregs; q++)
        initBlankState(quregs[q]);
        
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
}

void DerivCircuit::calcDerivEnergiesStateVec(qreal* eneryGrad, PauliHamil hamil, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (!circuit->isUnitary()) // throws
        throw QuESTException("", "The given circuit must be composed strictly of unitary gates "
            "(and ergo exclude gates like Matr[] and P[]) in order to return a valid real "
            "observable gradient. For non-unitary circuits, use ApplyCircuitDerivs[]."); // throws
    
    if (numWorkQuregs < 3)
        throw QuESTException("", "An internal error occured. Fewer than three working registers were "
            "passed to DerivCircuit::calcDerivEnergiesStateVec, despite prior validation."); // throws
        
    Qureg workLambda = workQuregs[0];
    Qureg workPhi    = workQuregs[1];
    Qureg workMu     = workQuregs[2];
    
    // prepare |lambda> = H circuit |init> and |phi> = circuit |init>
    cloneQureg(workLambda, initQureg);
    circuit->applyTo(workLambda); // throws
    cloneQureg(workPhi, workLambda);
    applyPauliHamil(workPhi, hamil, workLambda); // throws

    // clear energies
    for (int i=0; i<numVars; i++)
        eneryGrad[i] = 0;

    for (int t=numTerms-1; t>=0; t--) {
        
        DerivTerm derivTerm = terms[t];
        int gateInd = derivTerm.getGateInd();
        int varInd = derivTerm.getVarInd();
        
        // remove all gates >= gateInd (not removed by previous iteration) from workPhi
        circuit->applyInverseSubTo(workPhi, gateInd, 
            (t<numTerms-1)? terms[t+1].getGateInd() : circuit->getNumGates()); // throws

        // add all (daggered) gates > gateInd (not added by previous iteration) to workLambda
        circuit->applyDaggerSubTo(workLambda, gateInd + 1,
            (t<numTerms-1)? terms[t+1].getGateInd() + 1 : circuit->getNumGates()); // throws
            
        cloneQureg(workMu, workPhi);
        derivTerm.applyTo(workMu); // throws
        eneryGrad[varInd] += 2 * calcInnerProduct(workLambda, workMu).real;        
    }
}

void DerivCircuit::calcDerivEnergiesDenseHamil(qreal* energyGrad, Qureg hamilQureg, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (!circuit->isInvertible()) // throws
        throw QuESTException("", "The circuit must only contain invertible operators, and hence cannot "
            "contain measurements or projections. It is otherwise possible a general operator (like "
            "U or Kraus) was non-invertible for its particular parameter values. Please instead use ApplyCircuitDerivs[].");
        
    if (!circuit->isTracePreserving()) // throws
        throw QuESTException("", "The circuit must be trace-preserving and hence cannot contain operators "
            "like Fac[] and Matr[]. Please instead use KrausNonTP which tolerates numerical non-CPTP.");

    if (!extension_isHermitian(initQureg))
        throw QuESTException("", "The initial density matrix state must be Hermitian."); // throws
    
     if (numWorkQuregs < 3)
         throw QuESTException("", "An internal error occured. Fewer than three working registers were "
             "passed to DerivCircuit::calcDerivEnergiesDenseHamil, despite prior validation."); // throws

    Qureg leftQureg  = workQuregs[0];
    Qureg rightQureg = workQuregs[1];
    Qureg derivQureg = workQuregs[2];
    
    // clear energies
    for (int i=0; i<numVars; i++)
        energyGrad[i] = 0;
        
    // an optimisation to avoid applying then undoing gates beyond the final deriv
    int lastDerGateInd = terms[numTerms-1].getGateInd();

    // prepare ||right>> = Phi_(lastDerGateInd-1) ... Phi_0 ||in>>
    cloneQureg(rightQureg, initQureg);
    circuit->applySubTo(rightQureg, 0, lastDerGateInd);
    int indOfLastGateOnRight = lastDerGateInd - 1;

    // prepare <<left|| = <<hamil|| Phi_N ... Phi_(lastDerGateInd+1)
    cloneQureg(leftQureg, hamilQureg);
    circuit->applyDaggerSubTo(leftQureg, lastDerGateInd+1, circuit->getNumGates());
    int indOfLastGateOnLeft = lastDerGateInd+1;
        
    for (int t=numTerms-1; t>=0; t--) {
        
        DerivTerm derivTerm = terms[t];
        int gateInd = derivTerm.getGateInd();
        int varInd = derivTerm.getVarInd();
        
        // set ||right>> = Phi_(gateInd-1) ... Phi_0 ||in>>  (by gate 'undoing')
        circuit->applyInverseSubTo(rightQureg, gateInd, indOfLastGateOnRight+1);
        indOfLastGateOnRight = gateInd - 1;
        
        // set <<left|| = <<hamil|| Phi_N ... Phi_(gateInd+1) (by adding daggers)
        circuit->applyDaggerSubTo(leftQureg, gateInd+1, indOfLastGateOnLeft);
        indOfLastGateOnLeft = gateInd + 1;
        
        // set ||deriv>> = (d Phi_(gateInd)/dx) Phi_(gateInd-1) ... Phi_0 ||in>>
        cloneQureg(derivQureg, rightQureg);
        derivTerm.applyTo(derivQureg);
        
        // add <<left||qureg>> = Tr(hamil drho/dx)
        Complex energy = statevec_calcInnerProduct(leftQureg, derivQureg);
        if (local_isNonZero(energy.imag))
            throw QuESTException("", "The energy had a non-real component (" + local_qcompToStr(fromComplex(energy)) + 
                ") likely due to a non-trace-preserving operator (KrausNonTP) or an "
                "invalidly initialised pauli string Qureg (see ?SetQuregToPauliString).");
        
        energyGrad[varInd] += energy.real;
    }
}

void DerivCircuit::calcDerivEnergiesDensMatr(qreal* energyGrad, PauliHamil hamil, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (numWorkQuregs < 4)
        throw QuESTException("", "An internal error occured. Fewer than four working registers were "
            "passed to DerivCircuit::calcDerivEnergiesDensMatr (receiving PauliHamil, and hence "
            "intending to populate a register to a pauli string), despite prior validation."); // throws
    
    Qureg hamilQureg = workQuregs[3];
    setQuregToPauliHamil(hamilQureg, hamil);
    calcDerivEnergiesDenseHamil(energyGrad, hamilQureg, initQureg, workQuregs, 3);
}

void DerivCircuit::calcDerivEnergies(qreal* energyGrad, PauliHamil hamil, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (circuit->isPure() && !initQureg.isDensityMatrix) // throws
        calcDerivEnergiesStateVec(energyGrad, hamil, initQureg, workQuregs, numWorkQuregs); // throws
    else
        calcDerivEnergiesDensMatr(energyGrad, hamil, initQureg, workQuregs, numWorkQuregs); // throws
}

qmatrix DerivCircuit::calcMetricTensorStateVec(Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (!circuit->isInvertible()) // throws
        throw QuESTException("", "The circuit must only contain invertible operators, and hence cannot "
            "contain measurements or projections. Please instead use ApplyCircuitDerivs[]"); // throws
    
    if (numWorkQuregs < 4)
        throw QuESTException("", "An internal error occured. Fewer than four working registers were "
            "passed to DerivCircuit::calcMetricTensorStateVec, despite prior validation."); // throws
    
    Qureg quregDiag   = workQuregs[0];
    Qureg quregSuffix = workQuregs[1];
    Qureg quregPrefix = workQuregs[2];
    Qureg quregDeriv  = workQuregs[3];
    
    cloneQureg(quregDiag, initQureg);
    
    // clear tensor
    qmatrix tensor = local_getQmatrix(numVars);
    for (int i=0; i<numVars; i++)
        for (int j=0; j<numVars; j++)
            tensor[i][j] = 0;
            
    // separate array for the Berry connection terms
    qvector berries = local_getQvector(numVars);
    for (int i=0; i<numVars; i++)
        berries[i] = 0;
    
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
            
    return tensor;
}

qmatrix DerivCircuit::calcMetricTensorDensMatr(Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (!circuit->isInvertible()) // throws
        throw QuESTException("", "The circuit must only contain invertible operators, and hence cannot "
            "contain measurements or projections. It is otherwise possible a general operator (like "
            "U or Kraus) was non-invertible for its particular parameter values. Please instead use ApplyCircuitDerivs[].");
        
    if (!circuit->isTracePreserving()) // throws
        throw QuESTException("", "The circuit must be trace-preserving and hence cannot contain operators "
            "like Fac[] and Matr[]. Please instead use KrausNonTP which tolerates numerical non-CPTP.");

    if (!extension_isHermitian(initQureg))
        throw QuESTException("", "The initial density matrix state must be Hermitian."); // throws
    
     if (numWorkQuregs < 4)
         throw QuESTException("", "An internal error occured. Fewer than four working registers were "
             "passed to DerivCircuit::calcMetricTensorDensMatr, despite prior validation."); // throws
    
    Qureg diagQureg = workQuregs[0];
    Qureg leftQureg = workQuregs[1];
    Qureg rightQureg = workQuregs[2];
    Qureg derivQureg = workQuregs[3];
    
    qmatrix tensor = local_getQmatrix(numVars);
    for (int i=0; i<numVars; i++)
        for (int j=0; j<numVars; j++)
            tensor[i][j] = 0;
    
    // set ||diagQureg>> = Phi_(finalDerivInd -1) ... Phi_0 ||in>>
    cloneQureg(diagQureg, initQureg);
    int finDerGateInd = terms[numTerms-1].getGateInd();
    circuit->applySubTo(diagQureg, 0, finDerGateInd); // throws
    int indOfLastGateOnDiag = finDerGateInd-1;

    for (int t=numTerms-1; t>=0; t--) {
        
        DerivTerm leftDerivTerm = terms[t];
        int leftVarInd = leftDerivTerm.getVarInd();
        int leftGateInd = leftDerivTerm.getGateInd(); // ordered
            
        // create ||diag>> = Phi_(leftInd-1) ... Phi_0 ||in>>
        circuit->applyInverseSubTo(diagQureg, leftGateInd, indOfLastGateOnDiag+1);        
        indOfLastGateOnDiag = leftGateInd-1;

        // set <<left|| = <<in|| Phi_0^ ... deriv(Phi_leftInd)^ ... Phi_last^
        cloneQureg(leftQureg, diagQureg);
        leftDerivTerm.applyTo(leftQureg); // throws
        circuit->applySubTo(leftQureg, leftGateInd+1, circuit->getNumGates()); // throws
        
        // compute diagonal element, via Tr(dagger(left)right) = <left|right>
        cloneQureg(rightQureg, leftQureg);
        Complex prod = statevec_calcInnerProduct(leftQureg, rightQureg);
        tensor[leftVarInd][leftVarInd] += fromComplex(prod); // real
        
        // set ||right>> = Phi_(leftInd-1) ... Phi_0 ||in>>
        cloneQureg(rightQureg, diagQureg);
        
        int indOfLastGateOnLeft = circuit->getNumGates();
        int indOfLastGateOnRight = leftGateInd - 1;
        
        for (int s=t-1; s>=0; s--) {
            
            DerivTerm rightDerivTerm = terms[s];
            int rightGateInd = rightDerivTerm.getGateInd();
            int rightVarInd = rightDerivTerm.getVarInd();
    
            // set <<left|| = <<in|| Phi_0^ ... deriv(Phi_leftInd)^ ... Phi_last^ Phi_last ... Phi_(rightInd + 1)
            // (on the first iteration, this performs O(numVars) gates, but O(1) subsequently)
            circuit->applyDaggerSubTo(leftQureg, rightGateInd+1, indOfLastGateOnLeft); // throws
            indOfLastGateOnLeft = rightGateInd+1;
            
            // set ||right>> = Phi_(rightInd-1) ... Phi_0 ||in>> 
            circuit->applyInverseSubTo(rightQureg, rightGateInd, indOfLastGateOnRight+1); // throws
            indOfLastGateOnRight = rightGateInd - 1;
            
            // set ||deriv>> = deriv(Phi_rightInd) Phi_(rightInd-1) ... Phi_0 ||in>> 
            cloneQureg(derivQureg, rightQureg);
            rightDerivTerm.applyTo(derivQureg);
        
            // HS(left,right) = <<left||right>> = 
            // <<in|| Phi_0^ ... deriv(Phi_leftInd)^ ... Phi_last^ Phi_last ... deriv(Phi_rightInd) ... Phi_0 ||in>> 
            prod = statevec_calcInnerProduct(leftQureg, derivQureg);
            tensor[leftVarInd][rightVarInd] += fromComplex(prod);
            tensor[rightVarInd][leftVarInd] += conj(fromComplex(prod));
        }
    }

    return tensor;
}

qmatrix DerivCircuit::calcMetricTensor(Qureg initQureg, Qureg* workQuregs, int numWorkQuregs) {
    
    if (circuit->isPure() && !initQureg.isDensityMatrix) // throws
        return calcMetricTensorStateVec(initQureg, workQuregs, numWorkQuregs); // throws
    else
        return calcMetricTensorDensMatr(initQureg, workQuregs, numWorkQuregs); // throws
}

int DerivCircuit::getNumNeededWorkQuregsFor(std::string funcName, Qureg initQureg) {
    
    int circIsPure = circuit->isPure(); // throws
    
    if (funcName == "applyTo")
        return 1;
    
    if (funcName == "calcDerivEnergies") {
        
        if (circIsPure && !initQureg.isDensityMatrix)
            return 3;
        else
            return 4;
    }
    
    if (funcName == "calcDerivEnergiesDenseHamil") {
        
        return 3;
    }
    
    if (funcName == "calcMetricTensor") {
        
        if (circIsPure && !initQureg.isDensityMatrix)
            return 4;
        else
            return 4;
    }
    
    local_sendErrorAndFail("ApplyCircuitDerivs", "An irrevocable internal error occured (DerivCircuit::getNumNeededWorkQuregsFor() " 
        "was given an unrecognised funcName: " + funcName + ") and the QuESTlink process must crash.");
    throw QuESTException("", "An internal error occurred; the function named passed to getNumNeededWorkQuregsFor() was unrecognised."); // throws
}

void DerivCircuit::validateWorkQuregsFor(std::string methodName, int initQuregId, int* workQuregIds, int numWorkQuregs) {
    
    // no working registers is fine; they will be internally created
    if (numWorkQuregs == 0)
        return;
        
    int numNeeded = getNumNeededWorkQuregsFor(methodName, quregs[initQuregId]); // throws
    if (numWorkQuregs < numNeeded)
        throw QuESTException("", "Too few working registers were passed (" + std::to_string(numNeeded) + " are required)."); // throws
        
    for (int i=0; i<numWorkQuregs; i++) {
        
        int workspaceId = workQuregIds[i];
        local_throwExcepIfQuregNotCreated(workspaceId); // throws
        
        if (workspaceId == initQuregId) 
            throw QuESTException("", "The initial state qureg must not be included in the workspace quregs."); // throws
            
        int numQb = quregs[initQuregId].numQubitsRepresented;
        int isDens = quregs[initQuregId].isDensityMatrix;
        
        if (quregs[workspaceId].numQubitsRepresented != numQb || quregs[workspaceId].isDensityMatrix != isDens)
            throw QuESTException("", "Workspace quregs must have the same type and dimension as the initial state qureg."); // throws
    }        
}

DerivCircuit::~DerivCircuit() {
    
    freeMMA();
    delete circuit;
    delete[] terms;
}



/*
 * interfacing
 */

void internal_applyCircuitDerivs(int initQuregId, int workspaceId) {
    const std::string apiFuncName = "ApplyCircuitDerivs";
    
    // get qureg ids (one for each var)
    int* quregIds;
    int numQuregs;
    WSGetInteger32List(stdlink, &quregIds, &numQuregs);
    
    // load the circuit and derivative descriptions (local so no need to explicitly delete)
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA();
    
    // validate quregs (must do so after loading derivCircuit from MMA so those packets are flushed)
    try {
        // validate initial state
        local_throwExcepIfQuregNotCreated(initQuregId); // throws
        
        // validate workspace (if pre-allocated)
        if (workspaceId != -1) {
            int workspaces[] = {workspaceId}; // hacky integration into validateWorkQuregsFor(), due to design indecision
            derivCirc.validateWorkQuregsFor("applyTo", initQuregId, workspaces, 1);
        }
            
        // validate derivative registers (no check of uniqueness)
        int numQb = quregs[initQuregId].numQubitsRepresented;
        int isDens = quregs[initQuregId].isDensityMatrix;
            
        for (int q=0; q<numQuregs; q++) {
            local_throwExcepIfQuregNotCreated(quregIds[q]); // throws
            if (quregIds[q] == initQuregId)
                throw QuESTException("", "The initial state qureg must not also be one of the derivative quregs."); // throws
            if (quregIds[q] == workspaceId)
                throw QuESTException("", "The workspace qureg must not be one of the derivative quregs."); // throws
            if (quregs[quregIds[q]].numQubitsRepresented != numQb)
                throw QuESTException("", "All derivative quregs must have the same dimension as the initial state."); // throws
            if (quregs[quregIds[q]].isDensityMatrix != isDens)
                throw QuESTException("", "Quregs must be all state-vectors or all density-matrices"); // throws
        }

        if (isDens && !extension_isHermitian(quregs[initQuregId]))
            throw QuESTException("", "The initial density matrix state must be Hermitian."); // throws
        
    } catch (QuESTException& err) {
        local_sendErrorAndFail(apiFuncName, err.message);
        WSReleaseInteger32List(stdlink, quregIds, numQuregs);
        return;
    }
    
    Qureg initQureg = quregs[initQuregId];
    
    // optionally create workspace
    Qureg workspace;
    if (workspaceId == -1)
        workspace = createCloneQureg(initQureg, env);
    else
        workspace = quregs[workspaceId];
    
    // unpack deriv quregs
    Qureg* derivQuregs = (Qureg*) malloc(numQuregs * sizeof *derivQuregs);
    for (int q=0; q<numQuregs; q++)
        derivQuregs[q] = quregs[quregIds[q]];
        
    // attempt to compute derivatives and return initQuregId
    try {
        derivCirc.applyTo(derivQuregs, numQuregs, initQureg, workspace); // throws
        WSPutInteger(stdlink, initQuregId);
        
    } catch (QuESTException& err) {
        local_sendErrorAndFailOrAbortFromExcep(apiFuncName, err.thrower,  err.message);
    } 
    
    // clean-up even in event of error 
    free(derivQuregs);
    WSReleaseInteger32List(stdlink, quregIds, numQuregs);
    if (workspaceId == -1)
        destroyQureg(workspace, env);
}

void internal_calcExpecPauliStringDerivs(int initQuregId) {
    const std::string apiFuncName = "CalcExpecPauliStringDerivs";
    
    // load the any-length workspace list from MMA
    int* workQuregIds;
    int numPassedWorkQuregs;
    WSGetInteger32List(stdlink, &workQuregIds, &numPassedWorkQuregs);
    
    // load the circuit and deriv spec from MMA
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA(); // local, so desconstructor automatic
    
    // load Hamiltonian from MMA (and also validate initQuregId)
    PauliHamil hamil;
    try {
        hamil = local_loadPauliHamilForQuregFromMMA(initQuregId); // throws
        
    } catch (QuESTException& err) {
        local_sendErrorAndFail(apiFuncName, err.message);
        WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
        return;
    }
    
    // validate registers 
    try {
        if (numPassedWorkQuregs > 0)
            derivCirc.validateWorkQuregsFor("calcDerivEnergies", initQuregId, workQuregIds, numPassedWorkQuregs); // throws
            
    } catch (QuESTException& err) {
        local_sendErrorAndFail(apiFuncName, err.message);
        WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
        local_freePauliHamil(hamil);
        return;
    }
    
    Qureg initQureg = quregs[initQuregId];
    
    // optionally create work registers
    int numNeededWorkQuregs = derivCirc.getNumNeededWorkQuregsFor("calcDerivEnergies", initQureg);
    Qureg* workQuregs = (Qureg*) malloc(numNeededWorkQuregs * sizeof *workQuregs);
    for (int i=0; i<numNeededWorkQuregs; i++)
        if (numPassedWorkQuregs == 0)
            workQuregs[i] = createCloneQureg(initQureg, env);
        else
            workQuregs[i] = quregs[workQuregIds[i]];
            
    // prepare energy grad vector (malloc onto heap to avoid stack size limits)
    int numDerivs = derivCirc.getNumVars();
    qreal* energyGrad = (qreal*) malloc(numDerivs * sizeof *energyGrad);
    
    // attempt to compute and return energy grad
    try {    
        derivCirc.calcDerivEnergies(energyGrad, hamil, initQureg, workQuregs, numNeededWorkQuregs); // throws
        WSPutReal64List(stdlink, energyGrad, numDerivs);
        
    } catch (QuESTException& err) {
        local_sendErrorAndFailOrAbortFromExcep(apiFuncName, err.thrower,  err.message);
    }

    // clean-up even despite errors
    if (numPassedWorkQuregs == 0)
        for (int i=0; i<numNeededWorkQuregs; i++)
            destroyQureg(workQuregs[i], env);
    free(workQuregs);
    free(energyGrad);
    local_freePauliHamil(hamil);
    WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
}

void internal_calcExpecPauliStringDerivsDenseHamil(int initQuregId, int hamilQuregId) {
    const std::string apiFuncName = "CalcExpecPauliStringDerivs";
    
    // load the any-length workspace list from MMA
    int* workQuregIds;
    int numPassedWorkQuregs;
    WSGetInteger32List(stdlink, &workQuregIds, &numPassedWorkQuregs);
    
    // load the circuit and deriv spec from MMA
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA(); // local, so desconstructor automatic
    
    Qureg initQureg, hamilQureg;
    
    // validate registers 
    try {
        local_throwExcepIfQuregNotCreated(initQuregId); // throws
        local_throwExcepIfQuregNotCreated(hamilQuregId); // throws
        
        initQureg = quregs[initQuregId];
        hamilQureg = quregs[hamilQuregId];
        
        if ( ! (initQureg.isDensityMatrix && hamilQureg.isDensityMatrix) )
            throw QuESTException("", "When passing a Hamiltonian encoded into a Qureg, "
                "both that Qureg and the initial state must be density matrices.");
                
        if ( initQureg.numQubitsRepresented != hamilQureg.numQubitsRepresented )
            throw QuESTException("", "The initial state Qureg and the Hamiltonian encoded into a "
                "qureg must be of equal dimensions.");
        
        if (numPassedWorkQuregs > 0)
            derivCirc.validateWorkQuregsFor("calcDerivEnergiesDenseHamil", initQuregId, workQuregIds, numPassedWorkQuregs); // throws
        
    } catch (QuESTException& err) {
        local_sendErrorAndFail(apiFuncName, err.message);
        WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
        return;
    }
    
    // optionally create work registers
    int numNeededWorkQuregs = derivCirc.getNumNeededWorkQuregsFor("calcDerivEnergiesDenseHamil", initQureg);
    Qureg* workQuregs = (Qureg*) malloc(numNeededWorkQuregs * sizeof *workQuregs);
    for (int i=0; i<numNeededWorkQuregs; i++)
        if (numPassedWorkQuregs == 0)
            workQuregs[i] = createCloneQureg(initQureg, env);
        else
            workQuregs[i] = quregs[workQuregIds[i]];
    
    // prepare energy grad vector (malloc onto heap to avoid stack size limits)
    int numDerivs = derivCirc.getNumVars();
    qreal* energyGrad = (qreal*) malloc(numDerivs * sizeof *energyGrad);
    
    // attempt to compute and return energy grad
    try {    
        derivCirc.calcDerivEnergiesDenseHamil(energyGrad, hamilQureg, initQureg, workQuregs, numNeededWorkQuregs); // throws

        WSPutReal64List(stdlink, energyGrad, numDerivs);
        
    } catch (QuESTException& err) {
        local_sendErrorAndFailOrAbortFromExcep(apiFuncName, err.thrower,  err.message);
    }

    // clean-up even despite errors
    if (numPassedWorkQuregs == 0)
        for (int i=0; i<numNeededWorkQuregs; i++)
            destroyQureg(workQuregs[i], env);
    free(workQuregs);
    free(energyGrad);
    WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
}

void internal_calcMetricTensor(int initQuregId) {
    const std::string apiFuncName = "CalcMetricTensor";
    
    // load the any-length workspace list from MMA
    int* workQuregIds;
    int numPassedWorkQuregs;
    WSGetInteger32List(stdlink, &workQuregIds, &numPassedWorkQuregs);
    
    // load the circuit and deriv spec from MMA
    DerivCircuit derivCirc;
    derivCirc.loadFromMMA(); // local, so desconstructor automatic
    
    // validate registers 
    try {
        local_throwExcepIfQuregNotCreated(initQuregId); // throws
        if (numPassedWorkQuregs > 0)
            derivCirc.validateWorkQuregsFor("calcMetricTensor", initQuregId, workQuregIds, numPassedWorkQuregs); // throws
            
    } catch (QuESTException& err) {
        local_sendErrorAndFail(apiFuncName, err.message);
        WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
        return;
    }
    
    Qureg initQureg = quregs[initQuregId];
    
    // optionally create work registers
    int numNeededWorkQuregs = derivCirc.getNumNeededWorkQuregsFor("calcMetricTensor", initQureg);
    Qureg* workQuregs = (Qureg*) malloc(numNeededWorkQuregs * sizeof *workQuregs);
    for (int i=0; i<numNeededWorkQuregs; i++)
        if (numPassedWorkQuregs == 0)
            workQuregs[i] = createCloneQureg(initQureg, env);
        else
            workQuregs[i] = quregs[workQuregIds[i]];
        
    // attempt to compute and return tensor 
    try {
        qmatrix tensor = derivCirc.calcMetricTensor(initQureg, workQuregs, numNeededWorkQuregs); // throws
        local_sendMatrixToMMA(tensor);
        
    } catch (QuESTException& err) {
        local_sendErrorAndFailOrAbortFromExcep(apiFuncName, err.thrower,  err.message);
    }
    
    // clean-up, even if error
    if (numPassedWorkQuregs == 0)
        for (int i=0; i<numNeededWorkQuregs; i++)
            destroyQureg(workQuregs[i], env);
    free(workQuregs);
    WSReleaseInteger32List(stdlink, workQuregIds, numPassedWorkQuregs);
}
