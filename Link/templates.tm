:: @file
:: Wraps C++ functions with Mathematica-callable symbols.
::
:: Some wrappers herein are directly callable by users without any surrounding 
:: definition in QuESTlink.m, and so much define ::usage and ::error tags 
:: (via :Evaluate). These use the QuEST` context. 
:: Other wrappers are only internally called by QuESTlink.m, and use the 
:: QuEST`Private` context, and an 'Internal' suffix.
::
:: ReturnType is always 'Manual' so that in the event of a user-input error, 
:: errors can be propagated back to the Mathematica frontend.
::
:: @author Tyson Jones 



:Begin:
:Function:       wrapper_createQureg
:Pattern:        QuEST`CreateQureg[numQubits_Integer]
:Arguments:      { numQubits }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CreateQureg::usage = "CreateQureg[numQubits] returns the id of a newly created statevector.";
    QuEST`CreateQureg::error = "`1`";
    QuEST`CreateQureg[___] := QuEST`Private`invalidArgError[CreateQureg];

:Begin:
:Function:       wrapper_createDensityQureg
:Pattern:        QuEST`CreateDensityQureg[numQubits_Integer]
:Arguments:      { numQubits }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CreateDensityQureg::usage = "CreateDensityQureg[numQubits] returns the id of a newly created density matrix.";
    QuEST`CreateDensityQureg::error = "`1`";
    QuEST`CreateDensityQureg[___] := QuEST`Private`invalidArgError[CreateDensityQureg];

:Begin:
:Function:       wrapper_destroyQureg
:Pattern:        QuEST`Private`DestroyQuregInternal[id_Integer]
:Arguments:      { id }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`DestroyQuregInternal::usage = "DestroyQuregInternal[numQubits] frees the memory of the qureg associated with the given id."

:Begin:
:Function:       callable_createQuregs
:Pattern:        QuEST`CreateQuregs[numQubits_Integer, numQuregs_Integer]
:Arguments:      { numQubits, numQuregs }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CreateQuregs::usage = "CreateQuregs[numQubits, numQuregs] returns a list of ids of a newly created statevectors.";
    QuEST`CreateQuregs::error = "`1`";
    QuEST`CreateQuregs[___] := QuEST`Private`invalidArgError[CreateQuregs];

:Begin:
:Function:       callable_createDensityQuregs
:Pattern:        QuEST`CreateDensityQuregs[numQubits_Integer, numQuregs_Integer]
:Arguments:      { numQubits, numQuregs }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CreateDensityQuregs::usage = "CreateDensityQuregs[numQubits, numQuregs] returns a list of ids of a newly created density matrices.";
    QuEST`CreateDensityQuregs::error = "`1`";
    QuEST`CreateDensityQuregs[___] := QuEST`Private`invalidArgError[CreateDensityQuregs];




:Begin:
:Function:       callable_startRecordingQASM
:Pattern:        QuEST`StartRecordingQASM[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`StartRecordingQASM::usage = "StartRecordingQASM[qureg] initiates QASM recording of subsequent gates applied to qureg.";
    QuEST`StartRecordingQASM::error = "`1`";
    QuEST`StartRecordingQASM[___] := QuEST`Private`invalidArgError[StartRecordingQASM];

:Begin:
:Function:       callable_stopRecordingQASM
:Pattern:        QuEST`StopRecordingQASM[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`StopRecordingQASM::usage = "StopRecordingQASM[qureg] stops QASM recording of subsequent gates to qureg. This does not affect the QASM recorded so far.";
    QuEST`StopRecordingQASM::error = "`1`";
    QuEST`StopRecordingQASM[___] := QuEST`Private`invalidArgError[StopRecordingQASM];

:Begin:
:Function:       callable_clearRecordedQASM
:Pattern:        QuEST`ClearRecordedQASM[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`ClearRecordedQASM::usage = "ClearRecordedQASM[qureg] clears all QASM recorded for the given qureg so far, but does not stop recording of subsequent gates.";
    QuEST`ClearRecordedQASM::error = "`1`";
    QuEST`ClearRecordedQASM[___] := QuEST`Private`invalidArgError[ClearRecordedQASM];
    
:Begin:
:Function:       callable_getRecordedQASM
:Pattern:        QuEST`GetRecordedQASM[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`GetRecordedQASM::usage = "GetRecordedQASM[qureg] returns the QASM recorded on qureg as a string.";
    QuEST`GetRecordedQASM::error = "`1`";
    QuEST`GetRecordedQASM[___] := QuEST`Private`invalidArgError[GetRecordedQASM];




:Begin:
:Function:       wrapper_initZeroState
:Pattern:        QuEST`InitZeroState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`InitZeroState::usage = "InitZeroState[qureg] sets the qureg to state |0> (and returns the qureg id).";
    QuEST`InitZeroState::error = "`1`";
    QuEST`InitZeroState[___] := QuEST`Private`invalidArgError[InitZeroState];

:Begin:
:Function:       wrapper_initPlusState
:Pattern:        QuEST`InitPlusState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer  }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`InitPlusState::usage = "InitPlusState[qureg] sets the qureg to state |+> (and returns the qureg id).";
    QuEST`InitPlusState::error = "`1`";
    QuEST`InitPlusState[___] := QuEST`Private`invalidArgError[InitPlusState];

:Begin:
:Function:       wrapper_initClassicalState
:Pattern:        QuEST`InitClassicalState[qureg_Integer, state_Integer]
:Arguments:      { qureg, state }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`InitClassicalState::usage = "InitClassicalState[qureg, ind] sets the qureg to basis state |ind> (and returns the qureg id).";
    QuEST`InitClassicalState::error = "`1`";
    QuEST`InitClassicalState[___] := QuEST`Private`invalidArgError[InitClassicalState];

:Begin:
:Function:       wrapper_initPureState
:Pattern:        QuEST`InitPureState[targetQureg_Integer, pureQureg_Integer]
:Arguments:      { targetQureg, pureQureg }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate:
    QuEST`InitPureState::usage = "InitPureState[targetQureg, pureQureg] puts targetQureg (statevec or density matrix) into the pureQureg (statevec) state (and returns the targetQureg id).";
    QuEST`InitPureState::error = "`1`";
    QuEST`InitPureState[___] := QuEST`Private`invalidArgError[InitPureState];

:Begin:
:Function:       wrapper_initStateFromAmps
:Pattern:        QuEST`InitStateFromAmps[qureg_Integer, reals_List, imags_List]
:Arguments:      { qureg, reals, imags }
:ArgumentTypes:  { Integer, RealList, RealList }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`InitStateFromAmps::usage = "InitStateFromAmps[qureg, reals, imags] initialises the given qureg to have the supplied amplitudes (and returns the qureg id).";
    QuEST`InitStateFromAmps::error = "`1`";
    QuEST`InitStateFromAmps[___] := QuEST`Private`invalidArgError[InitStateFromAmps];

:Begin:
:Function:       wrapper_cloneQureg
:Pattern:        QuEST`CloneQureg[target_Integer, source_Integer]
:Arguments:      { target, source }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate:
    QuEST`CloneQureg::usage = "CloneQureg[dest, source] sets dest to be a copy of source.";
    QuEST`CloneQureg::error = "`1`";
    QuEST`CloneQureg[___] := QuEST`Private`invalidArgError[CloneQureg];

:Begin:
:Function:       internal_getAmp
:Pattern:        QuEST`Private`GetAmpInternal[qureg_Integer, row_Integer, col_Integer]
:Arguments:      { qureg, row, col }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`GetAmpInternal::usage = "GetAmpInternal[qureg, row, col] returns complex amplitude with index [row] in a statevector qureg, or index [row][col] of a density matrix."

:Begin:
:Function:       internal_setAmp
:Pattern:        QuEST`Private`SetAmpInternal[qureg_Integer, ampRe_Real, ampIm_Real, row_Integer, col_Integer]
:Arguments:      { qureg, ampRe, ampIm, row, col}
:ArgumentTypes:  { Integer, Real, Real, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`SetAmpInternal::usage = "SetAmpInternal[qureg, ampRe, ampIm, row, col] modifies the amplitude with index [row] in a statevector qureg, or index [row][col] of a density matrix, to amplitude (ampRe + i*ampIm)."


:Begin:
:Function:       callable_isDensityMatrix
:Pattern:        QuEST`IsDensityMatrix[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`IsDensityMatrix::usage = "IsDensityMatrix[qureg] returns 0 or 1 to indicate whether qureg is a statevector or density matrix (respectively).";
    QuEST`IsDensityMatrix::error = "`1`";
    QuEST`IsDensityMatrix[___] := QuEST`Private`invalidArgError[IsDensityMatrix];


:Begin:
:Function:       internal_applyPhaseFunc
:Pattern:        QuEST`Private`ApplyPhaseFuncInternal[quregId_Integer, qubits_List, encoding_Integer, coeffs_List, exponents_List, overrideInds_List, overridePhases_List]
:Arguments:      { quregId, qubits, encoding, coeffs, exponents, overrideInds, overridePhases }
:ArgumentTypes:  { Integer, IntegerList, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`Private`ApplyPhaseFuncInternal::usage = "ApplyPhaseFuncInternal[quregId, qubits, encoding, coeffs, exponents, overrideInds, overridePhases] applies a diagonal unitary operator upon the qureg, with elements informed by the exponential-polynomial encoded in {coeffs}, {exponents}, applied to the state index informed by {qubits}.";

:Begin:
:Function:       internal_applyMultiVarPhaseFunc
:Pattern:        QuEST`Private`ApplyMultiVarPhaseFuncInternal[quregId_Integer, qubits_List, numQubitsPerReg_List, encoding_Integer, coeffs_List, exponents_List, numTermsPerReg_List, overrideInds_List, overridePhases_List]
:Arguments:      { quregId, qubits, numQubitsPerReg, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`Private`ApplyMultiVarPhaseFuncInternal::usage = "ApplyMultiVarPhaseFuncInternal[quregId, qubits, numQubitsPerReg, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases] applies a diagonal unitary operator upon the qureg, with elements informed by the multi-variable exponential-polynomial encoded in {coeffs}, {exponents}, applied to the state indices informed by {qubits}.";

:Begin:
:Function:       internal_applyNamedPhaseFunc
:Pattern:        QuEST`Private`ApplyNamedPhaseFuncInternal[quregId_Integer, qubits_List, numQubitsPerReg_List, encoding_Integer, phaseFuncName_Integer, overrideInds_List, overridePhases_List]
:Arguments:      { quregId, qubits, numQubitsPerReg, encoding, phaseFuncName, overrideInds, overridePhases }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`Private`ApplyNamedPhaseFuncInternal::usage = "ApplyNamedPhaseFuncInternal[quregId, qubits, numQubitsPerReg, encoding, phaseFuncName, overrideInds, overridePhases] applies a diagonal unitary operator upon the qureg, with elements informed by the multi-variable function implied by phaseFuncName, applied to the state indices informed by {qubits}.";

:Begin:
:Function:       internal_applyParamNamedPhaseFunc
:Pattern:        QuEST`Private`ApplyParamNamedPhaseFuncInternal[quregId_Integer, qubits_List, numQubitsPerReg_List, encoding_Integer, phaseFuncName_Integer, params_List, overrideInds_List, overridePhases_List]
:Arguments:      { quregId, qubits, numQubitsPerReg, encoding, phaseFuncName, params, overrideInds, overridePhases }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`Private`ApplyParamNamedPhaseFuncInternal::usage = "ApplyParamNamedPhaseFuncInternal[quregId, qubits, numQubitsPerReg, encoding, phaseFuncName, params, overrideInds, overridePhases] applies a diagonal unitary operator upon the qureg, with elements informed by the multi-variable function implied by phaseFuncName, applied to the state indices informed by {qubits}.";


:Begin:
:Function:       wrapper_mixDepolarising
:Pattern:        QuEST`MixDepolarising[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`MixDepolarising::usage = "MixDepolarising[qureg, qubit, prob] adds depolarising noise to density matrix qureg.";
    QuEST`MixDepolarising::error = "`1`";
    QuEST`MixDepolarising[___] := QuEST`Private`invalidArgError[MixDepolarising];

:Begin:
:Function:       wrapper_mixTwoQubitDepolarising
:Pattern:        QuEST`MixTwoQubitDepolarising[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Integer, Integer, Integer, Real }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`MixTwoQubitDepolarising::usage = "MixTwoQubitDepolarising[qureg, qb1, qb2 prob] adds depolarising noise to density matrix qureg.";
    QuEST`MixTwoQubitDepolarising::error = "`1`";
    QuEST`MixTwoQubitDepolarising[___] := QuEST`Private`invalidArgError[MixTwoQubitDepolarising];

:Begin:
:Function:       wrapper_mixDephasing
:Pattern:        QuEST`MixDephasing[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`MixDephasing::usage = "MixDephasing[qureg, qubit, prob] adds dephasing noise to density matrix qureg.";
    QuEST`MixDephasing::error = "`1`";
    QuEST`MixDephasing[___] := QuEST`Private`invalidArgError[MixDephasing];

:Begin:
:Function:       wrapper_mixTwoQubitDephasing
:Pattern:        QuEST`MixTwoQubitDephasing[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Integer, Integer, Integer, Real }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`MixTwoQubitDephasing::usage = "MixTwoQubitDephasing[qureg, qb1, qb2 prob] adds dephasing noise to density matrix qureg.";
    QuEST`MixTwoQubitDephasing::error = "`1`";
    QuEST`MixTwoQubitDephasing[___] := QuEST`Private`invalidArgError[MixTwoQubitDephasing];

:Begin:
:Function:       wrapper_mixDamping
:Pattern:        QuEST`MixDamping[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`MixDamping::usage = "MixDamping[qureg, qubit, prob] applies amplitude damping with the given decay probability to density matrix qureg.";
    QuEST`MixDamping::error = "`1`";
    QuEST`MixDamping[___] := QuEST`Private`invalidArgError[MixDamping];




:Begin:
:Function:       wrapper_calcProbOfOutcome
:Pattern:        QuEST`CalcProbOfOutcome[qureg_Integer, qb_Integer, outcome_Integer]
:Arguments:      { qureg, qb, outcome }
:ArgumentTypes:  { Integer, Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcProbOfOutcome::usage = "CalcProbOfOutcome[qureg, qubit, outcome] returns the probability of measuring qubit in the given outcome.";
    QuEST`CalcProbOfOutcome::error = "`1`";
    QuEST`CalcProbOfOutcome[___] := QuEST`Private`invalidArgError[CalcProbOfOutcome];
    
:Begin:
:Function:       wrapper_calcProbOfAllOutcomes
:Pattern:        QuEST`CalcProbOfAllOutcomes[qureg_Integer, qubits_List]
:Arguments:      { qureg, qubits }
:ArgumentTypes:  { Integer, IntegerList }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcProbOfAllOutcomes::usage = "CalcProbOfAllOutcomes[qureg, qubits] returns the probabilities of every classical substate of the given list of qubits. The probabilities are ordered by their corresponding classical value (increasing), assuming qubits is given least to most significant.";
    QuEST`CalcProbOfAllOutcomes::error = "`1`";
    QuEST`CalcProbOfAllOutcomes[___] := QuEST`Private`invalidArgError[CalcProbOfAllOutcomes];

:Begin:
:Function:       wrapper_calcFidelity
:Pattern:        QuEST`CalcFidelity[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcFidelity::usage = "CalcFidelity[qureg1, qureg2] returns the fidelity between the given states.";
    QuEST`CalcFidelity::error = "`1`";
    QuEST`CalcFidelity[___] := QuEST`Private`invalidArgError[CalcFidelity];

:Begin:
:Function:       wrapper_calcInnerProduct
:Pattern:        QuEST`CalcInnerProduct[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcInnerProduct::usage = "CalcInnerProduct[qureg1, qureg2] returns the complex inner product between the given states.";
    QuEST`CalcInnerProduct::error = "`1`";
    QuEST`CalcInnerProduct[___] := QuEST`Private`invalidArgError[CalcInnerProduct];

:Begin:
:Function:       wrapper_calcDensityInnerProduct
:Pattern:        QuEST`CalcDensityInnerProduct[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcDensityInnerProduct::usage = "CalcDensityInnerProduct[qureg1, qureg2] returns the real Hilbert schmidt scalar product between two given density matrices.";
    QuEST`CalcDensityInnerProduct::error = "`1`";
    QuEST`CalcDensityInnerProduct[___] := QuEST`Private`invalidArgError[CalcDensityInnerProduct];

:Begin:
:Function:       wrapper_calcPurity
:Pattern:        QuEST`CalcPurity[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcPurity::usage = "CalcPurity[qureg] returns the purity of the given density matrix.";
    QuEST`CalcPurity::error = "`1`";
    QuEST`CalcPurity[___] := QuEST`Private`invalidArgError[CalcPurity];

:Begin:
:Function:       wrapper_calcTotalProb
:Pattern:        QuEST`CalcTotalProb[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcTotalProb::usage = "CalcTotalProb[qureg] returns the total probability (normalisation) of the statevector (sum of abs-squared of amplitudes) or density matrix (trace), which should be 1.";
    QuEST`CalcTotalProb::error = "`1`";
    QuEST`CalcTotalProb[___] := QuEST`Private`invalidArgError[CalcTotalProb];

:Begin:
:Function:       wrapper_calcHilbertSchmidtDistance
:Pattern:        QuEST`CalcHilbertSchmidtDistance[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CalcHilbertSchmidtDistance::usage = "CalcHilbertSchmidtDistance[qureg1, qureg2] returns the Hilbert-Schmidt distance (Frobenius norm of the difference) between the given density matrices.";
    QuEST`CalcHilbertSchmidtDistance::error = "`1`";
    QuEST`CalcHilbertSchmidtDistance[___] := QuEST`Private`invalidArgError[CalcHilbertSchmidtDistance];

:Begin:
:Function:       internal_calcQuregDerivs
:Pattern:        QuEST`Private`CalcQuregDerivsInternal[initStateId_Integer, quregIds_List, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, numTargsPerOp_List, params_List, numParamsPerOp_List, derivOpInds_List, derivVarInds_List, derivParams_List, numDerivParamsPerDerivOp_List]
:Arguments:      { initStateId, quregIds, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, derivOpInds, derivVarInds, derivParams, numDerivParamsPerDerivOp }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcQuregDerivsInternal::usage = "CalcQuregDerivsInternal[initStateId, quregIds, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, derivOpInds, derivVarInds, derivParams, numDerivParamsPerDerivOp] accepts a circuit (complete with rotation angles) and a nominated set of gates (by indices), sets each qureg to be the result of applying the derivative of the circuit w.r.t the nominated gates, upon the initial state."

:Begin:
:Function:       internal_calcExpecPauliSumDerivs
:Pattern:        QuEST`Private`CalcExpecPauliSumDerivsInternal[initStateId_Integer, isPureCirc_Integer, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, numTargsPerOp_List, params_List, numParamsPerOp_List, derivOpInds_List, derivVarInds_List, derivParams_List, numDerivParamsPerDerivOp_List, termCoeffs_List, allPauliCodes_List, allPauliTargets_List, numPaulisPerTerm_List]
:Arguments:      { initStateId, isPureCirc, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, derivOpInds, derivVarInds, derivParams, numDerivParamsPerDerivOp, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm }
:ArgumentTypes:  { Integer, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcExpecPauliSumDerivsInternal::usage = "CalcExpecPauliSumDerivsInternal[initStateId, isPureCirc, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, derivOpInds, derivVarInds, derivParams, numDerivParamsPerDerivOp, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm] accepts a circuit (complete with rotation angles), a derivative specification, and a Hamiltonian, and returns the energy gradient."

:Begin:
:Function:       internal_calcGeometricTensor
:Pattern:        QuEST`Private`CalcGeometricTensorInternal[initStateId_Integer, isPureCirc_Integer, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, numTargsPerOp_List, params_List, numParamsPerOp_List, derivOpInds_List, derivVarInds_List, derivParams_List, numDerivParamsPerDerivOp_List]
:Arguments:      { initStateId, isPureCirc, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, derivOpInds, derivVarInds, derivParams, numDerivParamsPerDerivOp }
:ArgumentTypes:  { Integer, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcGeometricTensorInternal::usage = "CalcGeometricTensor[initStateId, isPureCirc, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp, derivOpInds, derivVarInds, derivParams, numDerivParamsPerDerivOp] accepts a circuit and derivative terms and returns the corresponding geometric tensor."

:Begin:
:Function:       internal_calcInnerProductsMatrix
:Pattern:        QuEST`Private`CalcInnerProductsMatrixInternal[quregIds_List]
:Arguments:      { quregIds }
:ArgumentTypes:  { IntegerList }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcInnerProductsMatrixInternal::usage = "CalcInnerProductsMatrixInternal[quregIds] returns seperate 1D (appending rows) lists for the real and imag components of the matrix with ith-jth element CalcInnerProduct[quregIds[i], quregIds[j]]."

:Begin:
:Function:       internal_calcInnerProductsVector
:Pattern:        QuEST`Private`CalcInnerProductsVectorInternal[braId_Integer, ketIds_List]
:Arguments:      { braId, ketIds }
:ArgumentTypes:  { Integer, IntegerList }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcInnerProductsVectorInternal::usage = "CalcInnerProductsVectorInternal[braId, ketIds] returns a vector with jth element CalcInnerProduct[braId, ketIds[j]]."

:Begin:
:Function:       internal_calcDensityInnerProductsMatrix
:Pattern:        QuEST`Private`CalcDensityInnerProductsMatrixInternal[quregIds_List]
:Arguments:      { quregIds }
:ArgumentTypes:  { IntegerList }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcDensityInnerProductsMatrixInternal::usage = "CalcDensityInnerProductsMatrixInternal[quregIds] returns a flattened real matrix (appending rows) with ith-jth element CalcDensityInnerProduct[quregIds[i], quregIds[j]]."

:Begin:
:Function:       internal_calcDensityInnerProductsVector
:Pattern:        QuEST`Private`CalcDensityInnerProductsVectorInternal[rhoId_Integer, omegaIds_List]
:Arguments:      { rhoId, omegaIds }
:ArgumentTypes:  { Integer, IntegerList }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcDensityInnerProductsVectorInternal::usage = "CalcDensityInnerProductsVectorInternal[rhoId, omegaIds] returns a vector with jth element CalcDensityInnerProduct[braId, ketIds[j]]."



:Begin:
:Function:       internal_applyCircuit
:Pattern:        QuEST`Private`ApplyCircuitInternal[qureg_Integer, storeBackup_Integer, showProgress_Integer, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, numTargsPerOp_List, params_List, numParamsPerOp_List]
:Arguments:      { qureg, storeBackup, showProgress, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp }
:ArgumentTypes:  { Integer, Integer, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`ApplyCircuitInternal::usage = "ApplyCircuitInternal[qureg, storeBackup, showProgress, opcodes, ctrls, numCtrlsPerOps, targs, numTargsPerOp, params, numParamsPerOps] applies a circuit (decomposed into codes) to the given qureg."

:Begin:
:Function:       internal_calcExpecPauliProd
:Pattern:        QuEST`Private`CalcExpecPauliProdInternal[qureg_Integer, workspace_Integer, paulis_List, targets_List]
:Arguments:      { qureg, workspace, paulis, targets }
:ArgumentTypes:  { Integer, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcExpecPauliProdInternal::usage = "CalcExpecPauliProdInternal[qureg, workspace, paulis, targets] returns the expected value of the qureg under the given pauli product. workspace must be a Qureg of equal dimensions to qureg."

:Begin:
:Function:       internal_calcExpecPauliSum
:Pattern:        QuEST`Private`CalcExpecPauliSumInternal[qureg_Integer, workspace_Integer, termCoeffs_List, allPauliCodes_List, allPauliTargets_List, numPaulisPerTerm_List]
:Arguments:      { qureg, workspace, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm }
:ArgumentTypes:  { Integer, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcExpecPauliSumInternal::usage = "CalcExpecPauliSumInternal[qureg, workspace, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm] returns the expected value of the qureg under the given sum of Pauli products, specified as flat lists. workspace must be a Qureg of equal dimensions to qureg."

:Begin:
:Function:       internal_applyPauliSum
:Pattern:        QuEST`Private`ApplyPauliSumInternal[inQureg_Integer, outQureg_Integer, termCoeffs_List, allPauliCodes_List, allPauliTargets_List, numPaulisPerTerm_List]
:Arguments:      { inQureg, outQureg, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm }
:ArgumentTypes:  { Integer, Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`ApplyPauliSumInternal::usage = "ApplyPauliSumInternal[inQureg, outQureg, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm] modifies outQureg under the given sum of Pauli products, specified as flat lists. inQureg and outQureg must have the same type and equal dimensions."

:Begin:
:Function:       internal_calcPauliSumMatrix
:Pattern:        QuEST`Private`CalcPauliSumMatrixInternal[numQubits_Integer, termCoeffs_List, allPauliCodes_List, allPauliTargets_List, numPaulisPerTerm_List]
:Arguments:      { numQubits, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`CalcPauliSumMatrixInternal::usage = "CalcPauliSumMatrixInternal[numQubits, termCoeffs, allPauliCodes, allPauliTargets, numPaulisPerTerm] returns the action of applying the given sum of Pauli products (specified as flat lists) to every basis state."

:Begin:
:Function:       internal_getQuregMatrix
:Pattern:        QuEST`Private`GetQuregMatrixInternal[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`GetQuregMatrixInternal::usage = "GetQuregMatrixInternal[qureg] returns the underlying statevector associated with the given qureg (flat, even for density matrices)."

:Begin:
:Function:       internal_setWeightedQureg
:Pattern:        QuEST`Private`SetWeightedQuregInternal[facRe1_Real,facIm1_Real,qureg1_Integer, facRe2_Real,facIm2_Real,qureg2_Integer, facReOut_Real,facImOut_Real,quregOut_Integer]
:Arguments:      { facRe1,facIm1,qureg1, facRe2,facIm2,qureg2, facReOut,facImOut,quregOut }
:ArgumentTypes:  { Real, Real, Integer, Real, Real, Integer, Real, Real, Integer }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`SetWeightedQuregInternal::usage = "SetWeightedQuregInternal[facRe1,facIm1,qureg1, facRe2,facIm2,qureg2, facReOut,facImOut,quregOut] modifies quregOut to become (fac1 qureg1 + fac2 qureg2 + facOut qurgeOut)."

:Begin:
:Function:       wrapper_collapseToOutcome
:Pattern:        QuEST`CollapseToOutcome[qureg_Integer, qubit_Integer, outcome_Integer]
:Arguments:      { qureg, qubit, outcome }
:ArgumentTypes:  { Integer, Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`CollapseToOutcome::usage = "CollapseToOutcome[qureg, qubit, outcome] forces the target qubit to collapse to the given outcome.";
    QuEST`CollapseToOutcome::error = "`1`";
    QuEST`CollapseToOutcome[___] := QuEST`Private`invalidArgError[CollapseToOutcome];
    


:Begin:
:Function:       wrapper_applyFullQFT
:Pattern:        QuEST`ApplyQFT[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Begin:
:Function:       wrapper_applyQFT
:Pattern:        QuEST`ApplyQFT[qureg_Integer, qubits_List]
:Arguments:      { qureg, qubits }
:ArgumentTypes:  { Integer, IntegerList }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`ApplyQFT::usage = "ApplyQFT[qureg] applies the quantum Fourier transform circuit to the entire register.\nApplyQFT[qureg, qubits] applies the quantum Fourier transform circuit to the given qubits, assuming least-significant first.";
    QuEST`ApplyQFT::error = "`1`";
    QuEST`ApplyQFT[___] := QuEST`Private`invalidArgError[ApplyQFT];



:Begin:
:Function:       callable_destroyAllQuregs
:Pattern:        QuEST`DestroyAllQuregs[]
:Arguments:      { }
:ArgumentTypes:  { }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`DestroyAllQuregs::usage = "DestroyAllQuregs[] destroys all remote quregs.";
    QuEST`DestroyAllQuregs::error = "`1`";
    QuEST`DestroyAllQuregs[___] := QuEST`Private`invalidArgError[DestroyAllQuregs];

:Begin:
:Function:       callable_getAllQuregs
:Pattern:        QuEST`GetAllQuregs[]
:Arguments:      { }
:ArgumentTypes:  { }
:ReturnType:     Manual
:End:
:Evaluate: 
    QuEST`GetAllQuregs::usage = "GetAllQuregs[] returns all active quregs.";
    QuEST`GetAllQuregs::error = "`1`";
    QuEST`GetAllQuregs[___] := QuEST`Private`invalidArgError[GetAllQuregs];