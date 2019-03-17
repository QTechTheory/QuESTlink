
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
:Evaluate: QuEST`InitZeroState::usage = "InitZeroState[qureg] returns a state in |0>."

:Begin:
:Function:       wrapper_initPlusState
:Pattern:        QuEST`InitPlusState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer  }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitPlusState::usage = "InitPlusState[qureg] returns a state in |+>."

:Begin:
:Function:       wrapper_initClassicalState
:Pattern:        QuEST`InitClassicalState[qureg_Integer, state_Integer]
:Arguments:      { qureg, state }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitClassicalState::usage = "InitClassicalState[qureg, ind] returns a state in basis state |ind>."

:Begin:
:Function:       wrapper_initPureState
:Pattern:        QuEST`InitPureState[targetQureg_Integer, pureQureg_Integer]
:Arguments:      { targetQureg, pureQureg }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitPureState::usage = "InitPureState[targetQureg, pureQureg] puts targetQureg (statevec or density matrix) into the pureQureg (statevec) state."

:Begin:
:Function:       wrapper_initStateFromAmps
:Pattern:        QuEST`InitStateFromAmps[qureg_Integer, reals_List, imags_List]
:Arguments:      { qureg, reals, imags }
:ArgumentTypes:  { Integer, RealList, RealList }
:ReturnType:     Integer
:End:
:Evaluate: QuEST`InitStateFromAmps::usage = "InitStateFromAmps[qureg, reals, imags] initialises the given qureg to have the supplied amplitudes."

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
:Function:       internal_applyCircuit
:Pattern:        QuEST`Private`ApplyCircuitInternal[qureg_Integer, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, params_List, numParamsPerOp_List]
:Arguments:      { qureg, opcodes, ctrls, numCtrlsPerOp, targs, params, numParamsPerOp }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`ApplyCircuitInner::usage = "ApplyCircuitInner[qureg, opcodes, ctrls, numCtrlsPerOps, targs, params, numParamsPerOps] applies a circuit (decomposed into codes) to the given qureg."

:Begin:
:Function:       internal_getStateVec
:Pattern:        QuEST`Private`GetStateVecInternal[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Manual
:End:
:Evaluate: QuEST`Private`GetStateVecInternal::usage = "GetStateVecInternal[qureg] returns the underlying statevector associated with the given qureg (flat, even for density matrices)."




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