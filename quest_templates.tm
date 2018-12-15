:Begin:
:Function:       applyCircuit
:Pattern:        ApplyCircuitInner[opcodes_List, ctrls_List, targs_List, params_List, qureg_List]
:Arguments:      { opcodes, ctrls, targs, params, qureg }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: ApplyCircuitInner::usage = "ApplyCircuitInner[opcodes, ctrls, targs, params, qureg] applies a circuit (decomposed into codes) to the given qureg."




:Begin:
:Function:       wrapper_calcProbOfOutcome
:Pattern:        CalcProbOfOutcome[qureg_List, qb_Integer, outcome_Integer]
:Arguments:      { qureg, qb, outcome }
:ArgumentTypes:  { Manual }
:ReturnType:     Real
:End:
:Evaluate: CalcProbOfOutcome::usage = "CalcProbOfOutcome[qureg, qubit, outcome] returns the probability of measuring qubit in the given outcome."

:Begin:
:Function:       wrapper_calcFidelity
:Pattern:        CalcFidelity[qureg1_List, qureg2_List]
:Arguments:      { qureg1, qureg2 }
:ArgumentTypes:  { Manual }
:ReturnType:     Real
:End:
:Evaluate: CalcFidelity::usage = "CalcFidelity[qureg1, qureg2] returns the fidelity between the given states."




:Begin:
:Function:       wrapper_initZeroState
:Pattern:        InitZeroState[qureg_List]
:Arguments:      { qureg }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: InitZeroState::usage = "InitZeroState[qureg] returns a state in |0>."

:Begin:
:Function:       wrapper_initPlusState
:Pattern:        InitPlusState[qureg_List]
:Arguments:      { qureg }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: InitPlusState::usage = "InitPlusState[qureg] returns a state in |+>."

:Begin:
:Function:       wrapper_initClassicalState
:Pattern:        InitClassicalState[qureg_List, state_Integer]
:Arguments:      { qureg, state }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: InitClassicalState::usage = "InitClassicalState[qureg, ind] returns a state in basis state |ind>."

:Begin:
:Function:       wrapper_initPureState
:Pattern:        InitPureState[targetQureg_List, pureQureg_List]
:Arguments:      { targetQureg, pureQureg }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: InitCPureState::usage = "InitPureState[targetQureg, pureQureg] puts targetQureg (statevec or density matrix) into the pureQureg (statevec) state"




:Begin:
:Function:       wrapper_applyOneQubitDepolariseError
:Pattern:        ApplyOneQubitDepolariseError[qureg_List, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: ApplyOneQubitDepolariseError::usage = "ApplyOneQubitDepolariseError[qureg, qubit, prob] adds depolarising noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyTwoQubitDepolariseError
:Pattern:        ApplyTwoQubitDepolariseError[qureg_List, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: ApplyTwoQubitDepolariseError::usage = "ApplyTwoQubitDepolariseError[qureg, qb1, qb2 prob] adds depolarising noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyOneQubitDephaseError
:Pattern:        ApplyOneQubitDephaseError[qureg_List, qb_Integer, prob_Real]
:Arguments:      { qureg, qb, prob }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: ApplyOneQubitDephaseError::usage = "ApplyOneQubitDephaseError[qureg, qubit, prob] adds dephasing noise to density matrix qureg."

:Begin:
:Function:       wrapper_applyTwoQubitDephaseError
:Pattern:        ApplyTwoQubitDephaseError[qureg_List, qb1_Integer, qb2_Integer, prob_Real]
:Arguments:      { qureg, qb1, qb2, prob }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: ApplyTwoQubitDephaseError::usage = "ApplyTwoQubitDephaseError[qureg, qb1, qb2 prob] adds dephasing noise to density matrix qureg."
