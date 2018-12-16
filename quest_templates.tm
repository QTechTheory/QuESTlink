
:Begin:
:Function:       wrapper_createQureg
:Pattern:        CreateQureg[numQubits_Integer]
:Arguments:      { numQubits }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: CreateQureg::usage = "CreateQureg[numQubits] returns the id of a newly created remote statevector."

:Begin:
:Function:       wrapper_createDensityQureg
:Pattern:        CreateDensityQureg[numQubits_Integer]
:Arguments:      { numQubits }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: CreateDensityQureg::usage = "CreateDensityQureg[numQubits] returns the id of a newly created remote density matrix."

:Begin:
:Function:       wrapper_destroyQureg
:Pattern:        DestroyQuregInner[id_Integer]
:Arguments:      { id }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: DestroyQureg::usage = "DestroyQureg[numQubits] frees the memory of the remote qureg associated with the given id."




:Begin:
:Function:       wrapper_initZeroState
:Pattern:        InitZeroState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer }
:ReturnType:     Integer
:End:
:Evaluate: InitZeroState::usage = "InitZeroState[qureg] returns a state in |0>."

:Begin:
:Function:       wrapper_initPlusState
:Pattern:        InitPlusState[qureg_Integer]
:Arguments:      { qureg }
:ArgumentTypes:  { Integer  }
:ReturnType:     Integer
:End:
:Evaluate: InitPlusState::usage = "InitPlusState[qureg] returns a state in |+>."

:Begin:
:Function:       wrapper_initClassicalState
:Pattern:        InitClassicalState[qureg_Integer, state_Integer]
:Arguments:      { qureg, state }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: InitClassicalState::usage = "InitClassicalState[qureg, ind] returns a state in basis state |ind>."

:Begin:
:Function:       wrapper_initPureState
:Pattern:        InitPureState[targetQureg_Integer, pureQureg_Integer]
:Arguments:      { targetQureg, pureQureg }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: InitCPureState::usage = "InitPureState[targetQureg, pureQureg] puts targetQureg (statevec or density matrix) into the pureQureg (statevec) state"




:Begin:
:Function:       wrapper_applyOneQubitDepolariseError
:Pattern:        ApplyOneQubitDepolariseError[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: ApplyOneQubitDepolariseError::usage = "ApplyOneQubitDepolariseError[qureg, qubit, prob] adds depolarising noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyTwoQubitDepolariseError
:Pattern:        ApplyTwoQubitDepolariseError[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Integer, Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: ApplyTwoQubitDepolariseError::usage = "ApplyTwoQubitDepolariseError[qureg, qb1, qb2 prob] adds depolarising noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyOneQubitDephaseError
:Pattern:        ApplyOneQubitDephaseError[qureg_Integer, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: ApplyOneQubitDephaseError::usage = "ApplyOneQubitDephaseError[qureg, qubit, prob] adds dephasing noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyTwoQubitDephaseError
:Pattern:        ApplyTwoQubitDephaseError[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Integer, Integer, Integer, Real }
:ReturnType:     Integer
:End:
:Evaluate: ApplyTwoQubitDephaseError::usage = "ApplyTwoQubitDephaseError[qureg, qb1, qb2 prob] adds dephasing noise to density matrix qureg."




:Begin:
:Function:       wrapper_calcProbOfOutcome
:Pattern:        CalcProbOfOutcome[qureg_Integer, qb_Integer, outcome_Integer]
:Arguments:      { qureg, qb, outcome }
:ArgumentTypes:  { Integer, Integer, Integer }
:ReturnType:     Real
:End:
:Evaluate: CalcProbOfOutcome::usage = "CalcProbOfOutcome[qureg, qubit, outcome] returns the probability of measuring qubit in the given outcome."

:Begin:
:Function:       wrapper_calcFidelity
:Pattern:        CalcFidelity[qureg1_Integer, qureg2_Integer]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Real
:End:
:Evaluate: CalcFidelity::usage = "CalcFidelity[qureg1, qureg2] returns the fidelity between the given states."




:Begin:
:Function:       applyCircuit
:Pattern:        ApplyCircuitInner[qureg_Integer, opcodes_List, ctrls_List, targs_List, params_List]
:Arguments:      { qureg, opcodes, ctrls, targs, params }
:ArgumentTypes:  { Integer, Manual }
:ReturnType:     Integer
:End:
:Evaluate: ApplyCircuitInner::usage = "ApplyCircuitInner[qureg, opcodes, ctrls, targs, params] applies a circuit (decomposed into codes) to the given qureg."
