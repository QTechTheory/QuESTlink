
:Begin:
:Function:       wrapper_createQureg
:Pattern:        QuEST`CreateQureg[numQubits_Integer]
:Arguments:      { numQubits }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`CreateQureg::usage = "CreateQureg[numQubits] returns the id of a newly created remote statevector."

:Begin:
:Function:       wrapper_createDensityQureg
:Pattern:        QuEST`CreateDensityQureg[numQubits_Integer]
:Arguments:      { numQubits }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`CreateDensityQureg::usage = "CreateDensityQureg[numQubits] returns the id of a newly created remote density matrix."

:Begin:
:Function:       wrapper_destroyQureg
:Pattern:        QuEST`Private`DestroyQuregInternal[id_Integer]
:Arguments:      { id }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`Private`DestroyQuregInternal::usage = "DestroyQuregInternal[numQubits] frees the memory of the remote qureg associated with the given id."




:Begin:
:Function:       wrapper_initZeroState
:Pattern:        QuEST`InitZeroState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitZeroState::usage = "InitZeroState[qureg] sets the qureg to state |0> (and returns the qureg id)."

:Begin:
:Function:       wrapper_initPlusState
:Pattern:        QuEST`InitPlusState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer  }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitPlusState::usage = "InitPlusState[qureg] sets the qureg to state |+> (and returns the qureg id)."

:Begin:
:Function:       wrapper_initClassicalState
:Pattern:        QuEST`InitClassicalState[qureg_Integer, state_Integer]
:Arguments:      { qureg, state }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitClassicalState::usage = "InitClassicalState[qureg, ind] sets the qureg to basis state |ind> (and returns the qureg id)."

:Begin:
:Function:       wrapper_initPureState
:Pattern:        QuEST`InitPureState[targetQureg_Integer, pureQureg_Integer]
:Arguments:      { targetQureg, pureQureg }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitPureState::usage = "InitPureState[targetQureg, pureQureg] puts targetQureg (statevec or density matrix) into the pureQureg (statevec) state (and returns the targetQureg id)."

:Begin:
:Function:       wrapper_initStateFromAmps
:Pattern:        QuEST`InitStateFromAmps[qureg_Integer, reals_List, imags_List]
:Arguments:      { qureg, reals, imags }
:ArgumentTypes:  { Integer, RealList, RealList }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitStateFromAmps::usage = "InitStateFromAmps[qureg, reals, imags] initialises the given qureg to have the supplied amplitudes (and returns the qureg id)."

:Begin:
:Function:       wrapper_cloneQureg
:Pattern:        QuEST`CloneQureg[target_Integer, source_Integer]
:Arguments:      { target, source }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`CloneQureg::usage = "CloneQureg[dest, source] sets dest to be a copy of source."




:Begin:
:Function:       wrapper_applyOneQubitDepolariseError
:Pattern:        QuEST`ApplyOneQubitDepolariseError[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`ApplyOneQubitDepolariseError::usage = "ApplyOneQubitDepolariseError[qureg, qubit, prob] adds depolarising noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyTwoQubitDepolariseError
:Pattern:        QuEST`ApplyTwoQubitDepolariseError[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Integer, Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`ApplyTwoQubitDepolariseError::usage = "ApplyTwoQubitDepolariseError[qureg, qb1, qb2 prob] adds depolarising noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyOneQubitDephaseError
:Pattern:        QuEST`ApplyOneQubitDephaseError[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`ApplyOneQubitDephaseError::usage = "ApplyOneQubitDephaseError[qureg, qubit, prob] adds dephasing noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyTwoQubitDephaseError
:Pattern:        QuEST`ApplyTwoQubitDephaseError[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Integer, Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`ApplyTwoQubitDephaseError::usage = "ApplyTwoQubitDephaseError[qureg, qb1, qb2 prob] adds dephasing noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyOneQubitDampingError
:Pattern:        QuEST`ApplyOneQubitDampingError[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`ApplyOneQubitDampingError::usage = "ApplyOneQubitDampingError[qureg, qubit, prob] applies amplitude damping with the given decay probability to density matrix qureg."




:Begin:
:Function:       wrapper_calcProbOfOutcome
:Pattern:        QuEST`CalcProbOfOutcome[qureg_Integer, qb_Integer, outcome_Integer]
:Arguments:      { qureg, qb, outcome }
:ArgumentTypes:  { Integer, Integer, Integer }
:ReturnType:     Real
:End:
:Evaluate: QuEST`CalcProbOfOutcome::usage = "CalcProbOfOutcome[qureg, qubit, outcome] returns the probability of measuring qubit in the given outcome."

:Begin:
:Function:       wrapper_calcFidelity
:Pattern:        QuEST`CalcFidelity[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Real
:End:
:Evaluate: QuEST`CalcFidelity::usage = "CalcFidelity[qureg1, qureg2] returns the fidelity between the given states."

:Begin:
:Function:       wrapper_calcInnerProduct
:Pattern:        QuEST`CalcInnerProduct[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`CalcInnerProduct::usage = "CalcInnerProduct[qureg1, qureg2] returns the complex inner product between the given states."

:Begin:
:Function:       wrapper_calcPurity
:Pattern:        QuEST`CalcPurity[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Real
:End:
:Evaluate: QuEST`CalcPurity::usage = "CalcPurity[qureg] returns the purity of the given density matrix."

:Begin:
:Function:       wrapper_calcTotalProb
:Pattern:        QuEST`CalcTotalProb[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Real
:End:
:Evaluate: QuEST`CalcTotalProb::usage = "CalcTotalProb[qureg] returns the total probability (normalisation) of the statevector (sum of abs-squared of amplitudes) or density matrix (trace), which should be 1."

:Begin:
:Function:       wrapper_calcHilbertSchmidtDistance
:Pattern:        QuEST`CalcHilbertSchmidtDistance[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Real
:End:
:Evaluate: QuEST`CalcHilbertSchmidtDistance::usage = "CalcHilbertSchmidtDistance[qureg1, qureg2] returns the Hilbert-Schmidt distance (Frobenius norm of the diference) between the given density matrices."




:Begin:
:Function:       internal_applyCircuit
:Pattern:        QuEST`Private`ApplyCircuitInternal[qureg_Integer, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, numTargsPerOp_List, params_List, numParamsPerOp_List]
:Arguments:      { qureg, opcodes, ctrls, numCtrlsPerOp, targs, numTargsPerOp, params, numParamsPerOp }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`ApplyCircuitInner::usage = "ApplyCircuitInner[qureg, opcodes, ctrls, numCtrlsPerOps, targs, numTargsPerOp, params, numParamsPerOps] applies a circuit (decomposed into codes) to the given qureg."

:Begin:
:Function:       internal_calcExpectedValue
:Pattern:        QuEST`Private`CalcExpectedValueInternal[qureg_Integer, paulis_List, targets_List]
:Arguments:      { qureg, paulis, targets }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Real
:End:
:Evaluate: QuEST`Private`CalcExpectedValueInternal::usage = "CalcExpectedValueInternal[qureg, paulis, targets] returns the expected value of the qureg under the given pauli product."

:Begin:
:Function:       internal_getStateVec
:Pattern:        QuEST`Private`GetStateVecInternal[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`GetStateVecInternal::usage = "GetStateVecInternal[qureg] returns the underlying statevector associated with the given qureg (flat, even for density matrices)."

:Begin:
:Function:       internal_addWeightedStates
:Pattern:        QuEST`Private`AddWeightedStatesInternal[facRe1_Real,facIm1_Real,qureg1_Integer, facRe2_Real,facIm2_Real,qureg2_Integer, facReOut_Real,facImOut_Real,quregOut_Integer]
:Arguments:      { facRe1,facIm1,qureg1, facRe2,facIm2,qureg2, facReOut,facImOut,quregOut }
:ArgumentTypes:  { Real, Real, Integer, Real, Real, Integer, Real, Real, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`Private`AddWeightedStatesInternal::usage = "AddWeightedStatesInternal[facRe1,facIm1,qureg1, facRe2,facIm2,qureg2, facReOut,facImOut,quregOut] modifies quregOut to become (fac1 qureg1 + fac2 qureg2 + facOut qurgeOut)."

:Begin:
:Function:       wrapper_collapseToOutcome
:Pattern:        QuEST`CollapseToOutcome[qureg_Integer, qubit_Integer, outcome_Integer]
:Arguments:      { qureg, qubit, outcome }
:ArgumentTypes:  { Integer, Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`CollapseToOutcome::usage = "CollapseToOutcome[qureg, qubit, outcome] forces the target qubit to collapse to the given outcome."




:Begin:
:Function:       callable_destroyAllQuregs
:Pattern:        QuEST`DestroyAllQuregs[]
:Arguments:      { }
:ArgumentTypes:  { }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`DestroyAllQuregs::usage = "DestroyAllQuregs[] destroys all remote quregs."

:Begin:
:Function:       callable_getAllQuregs
:Pattern:        QuEST`GetAllQuregs[]
:Arguments:      { }
:ArgumentTypes:  { }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`GetAllQuregs::usage = "GetAllQuregs[] returns all active quregs."