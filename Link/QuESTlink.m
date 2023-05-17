
(* @file
 * The front-end Mathematica interface between the user API and the backend C++ facilities.
 * Some functions herein merely wrap a core QuEST function, while others require 
 * complicated argument translation, and some invoke only Mathematica routines.
 * The C++ functions are wrapped to become callable Mathematica symbols via templates.tm
 *
 * @author Tyson Jones
 *)

BeginPackage["QuEST`"]    
    
    (* 
     * Note additional functions and their usage messages are fetched when CreateRemoteQuESTEnv is called.
     * The additional functions are provided through WSTP and use the QuEST` prefix to share namespace with this package.
     * This includes QuEST`Private` functions which are thus only called from within this package, which does not need to
     * explicitly use that package. 
     * Note also that public WSTP functions e.g. CloneQureg called from this package must be prefixed here with QuEST` so
     * that they are not incorrectly automatically given a QuEST`Private` prefix. This often will not trigger an error 
     * but will cause incorrect behaviour. When in doubt, give the full explicit name e.g. QuEST`CloneQureg or
     * QuEST`Private`ApplyCircuitInternal. 
     *)
     
     (*
      * public API 
      *)
    
    ApplyCircuit::usage = "ApplyCircuit[qureg, circuit] modifies qureg by applying the circuit. Returns any measurement outcomes and the probabilities encountered by projectors, ordered and grouped by the appearance of M and P in the circuit.
ApplyCircuit[inQureg, circuit, outQureg] leaves inQureg unchanged, but modifies outQureg to be the result of applying the circuit to inQureg.
Accepts optional arguments WithBackup and ShowProgress."
    ApplyCircuit::error = "`1`"
    
    ApplyCircuitDerivs::usage = "ApplyCircuitDerivs[inQureg, circuit, varVals, outQuregs] modifies outQuregs to be the result of applying the derivatives (with respect to variables in varVals) of the given symbolic circuit to inQureg (which remains unmodified).
    \[Bullet] varVals is a list {symbol -> value, ...} of all variables present in the circuit parameters.
    \[Bullet] outQuregs is a list of quregs to set to the respective derivative of circuit upon inQureg, according to the order of vars.
    \[Bullet] Variable repetition, multi-parameter gates, variable-dependent element-wise matrices, variable-dependent channels, and operators whose parameters are (numerically evaluable) functions of variables are all permitted within the circuit. In effect, every continuously-parameterised circuit or channel is permitted.
ApplyCircuitDerivs[inQureg, circuit, varVals, outQuregs, workQuregs] use the given persistent workspace quregs to avoid tediously creating and destroying any internal quregs, for a speedup. For convenience, any number of workspaces can be passed, but only the first is needed and used."
    ApplyCircuitDerivs::error = "`1`"
    
    CalcExpecPauliStringDerivs::usage = "CalcExpecPauliStringDerivs[inQureg, circuit, varVals, pauliString] returns the gradient vector of the pauliString expected values, as produced by the derivatives of the circuit (with respect to varVals, {var -> value}) acting upon the given initial state (inQureg).
CalcExpecPauliStringDerivs[inQureg, circuit, varVals, pauliQureg] accepts a Qureg pre-initialised as a pauli string via SetQuregToPauliString[] to speedup density-matrix simulation.
CalcExpecPauliStringDerivs[inQureg, circuit, varVals, pauliStringOrQureg, workQuregs] uses the given persistent workspaces (workQuregs) in lieu of creating them internally, and should be used for optimum performance. At most four workQuregs are needed.
    \[Bullet] Variable repetition, multi-parameter gates, variable-dependent element-wise matrices, variable-dependent channels, and operators whose parameters are (numerically evaluable) functions of variables are all permitted. 
    \[Bullet] All operators must be invertible, trace-preserving and deterministic, else an error is thrown. 
    \[Bullet] This function runs asymptotically faster than ApplyCircuitDerivs[] and requires only a fixed memory overhead."
    CalcExpecPauliStringDerivs::error = "`1`"
    
    CalcMetricTensor::usage = "CalcMetricTensor[inQureg, circuit, varVals] returns the natural gradient metric tensor, capturing the circuit derivatives (produced from initial state inQureg) with respect to varVals, specified with values {var -> value, ...}.
    CalcMetricTensor[inQureg, circuit, varVals, workQuregs] uses the given persistent workspace quregs (workQuregs) in lieu of creating them internally, and should be used for optimum performance. At most four workQuregs are needed.
    \[Bullet] For state-vectors and pure circuits, this returns the quantum geometric tensor, which relates to the Fubini-Study metric, the classical Fisher information matrix, and the variational imaginary-time Li tensor with Berry connections.
    \[Bullet] For density-matrices and noisy channels, this function returns the Hilbert-Schmidt derivative metric, which well approximates the quantum Fisher information matrix, though is a more experimentally relevant minimisation metric (https://arxiv.org/abs/1912.08660).
    \[Bullet] Variable repetition, multi-parameter gates, variable-dependent element-wise matrices, variable-dependent channels, and operators whose parameters are (numerically evaluable) functions of variables are all permitted. 
    \[Bullet] All operators must be invertible, trace-preserving and deterministic, else an error is thrown. 
    \[Bullet] This function runs asymptotically faster than ApplyCircuitDerivs[] and requires only a fixed memory overhead."
    CalcMetricTensor::error = "`1`"
    
    CalcInnerProducts::usage = "CalcInnerProducts[quregIds] returns a Hermitian matrix with i-th j-th element CalcInnerProduct[quregIds[i], quregIds[j]].
CalcInnerProducts[braId, ketIds] returns a complex vector with i-th element CalcInnerProduct[braId, ketIds[i]]."
    CalcInnerProducts::error = "`1`"

    CalcDensityInnerProducts::usage = "CalcDensityInnerProducts[quregIds] returns a Hermitian matrix with i-th j-th element CalcDensityInnerProduct[quregIds[i], quregIds[j]].
CalcDensityInnerProducts[rhoId, omegaIds] returns a vector with i-th element CalcDensityInnerProduct[rhoId, omegaIds[i]].
If all quregs are valid density matrices, the resulting tensors are real, though may have tiny non-zero imaginary components due to numerical imprecision.
For unnormalised density matrices, the tensors may contain complex scalars."
    CalcDensityInnerProducts::error = "`1`"
    
    Circuit::usage = "Circuit[gates] converts a product of gates into a left-to-right circuit, preserving order."
    Circuit::error = "`1`"
    
    Operator::usage = "Operator[gates] converts a product of gates into a right-to-left circuit."
    Operator::error = "`1`"

    CalcExpecPauliString::usage = "CalcExpecPauliString[qureg, pauliString, workspace] evaluates the expected value of a weighted sum of Pauli tensors, of a normalised qureg. workspace must be a qureg of equal dimensions to qureg. qureg is unchanged, and workspace is modified."
    CalcExpecPauliString::error = "`1`"

    ApplyPauliString::usage = "ApplyPauliString[inQureg, pauliString, outQureg] modifies outQureg to be the result of applying the weighted sum of Pauli tensors to inQureg."
    ApplyPauliString::error = "`1`"
    
    ApplyPhaseFunc::usage = "ApplyPhaseFunc[qureg, qubits, f[r], r] multiplies a phase factor e^(i f[r]) onto each amplitude in qureg, where r is substituted with the index of each basis state as informed by the list of qubits (ordered least to most significant), and optional argument BitEncoding.
\[Bullet] qubits is a list of which qubits to include in the determination of the index r for each basis state. For example, qubits={0,1,2} implies the canonical indexing of basis states in a 3-qubit register.
\[Bullet] f[r] must be an exponential polynomial of r, of the form sum_i a_j r^(p_j) where a_j and p_j can be any real number (including negative and fractional).
\[Bullet] f[r] must evaluate to a real number for every basis state index informed by qubits, unless overriden via optional argument PhaseOverrides.
ApplyPhaseFunc[qureg, {qubits, ...}, f[x,y,...], {x,y,...}] evaluates a multi-variable exponential-polynomial phase function, where each variable corresponds to a sub-register of qubits.
ApplyPhaseFunc[qureg, {qubits, ...}, FuncName] evaluates a specific named multi-variable function to determine the phase. These are:
    \[Bullet] \"Norm\" evaluates Sqrt[x^2 + y^2 + ...]
    \[Bullet] {\"InverseNorm\", div} evaluates 1/Sqrt[x^2 + y^2 + ...], replaced by div at divergence (when x=y=...=0).    
    \[Bullet] {\"ScaledNorm\", coeff} evaluates coeff*Sqrt[x^2 + y^2 + ...]
    \[Bullet] {\"ScaledInverseNorm\", coeff, div} evaluates coeff/Sqrt[x^2 + y^2 + ...], replaced by div at divergence (when x=y=...=0). 
    \[Bullet] {\"ScaledInverseShiftedNorm\", coeff, div, \[CapitalDelta]x, \[CapitalDelta]y, ...} evaluates coeff/Sqrt[(x-\[CapitalDelta]x)^2 + (y-\[CapitalDelta]y)^2 + ...], replaced by div at numerical divergence (when the denominator is within machine epsilon to zero). 
    \[Bullet] \"Product\" evaluates x*y*...
    \[Bullet] {\"InverseProduct\", div} evaluates 1/(x*y*...), replaced by div at divergence (when any of x, y, ... = 0).
    \[Bullet] {\"ScaledProduct\", coeff} evaluates coeff*x*y* ...
    \[Bullet] {\"ScaledInverseProduct\", coeff, div} evaluates coeff/(x*y* ...),, replaced by div at divergence (when any of x, y, ... = 0).
    \[Bullet] \"Distance\" evaluates Sqrt[(x1-x2)^2 + (y1-y2)^2 + ...], where sub-registers in {qubits} are assumed to be in order of {x1, x2, y1, y2, ...}
    \[Bullet] {\"InverseDistance\", div} evaluates 1/Sqrt[(x1-x2)^2 + (y1-y2)^2 + ...], replaced by div at divergence (when x1=x2, y1=y2, ...).   
    \[Bullet] {\"ScaledDistance\", coeff} evaluates coeff*Sqrt[(x1-x2)^2 + (y1-y2)^2 + ...]
    \[Bullet] {\"ScaledInverseDistance\", coeff, div} evaluates coeff/Sqrt[(x1-x2)^2 + (y1-y2)^2 + ...], replaced by div at divergence (when x1=x2, y1=y2, ...). 
    \[Bullet] {\"ScaledInverseShiftedDistance\", coeff, div, \[CapitalDelta]x, \[CapitalDelta]y, ...} evaluates coeff/Sqrt[(x1-x2-\[CapitalDelta]x)^2 + (y1-y2-\[CapitalDelta]y)^2 + ...], replaced by div at numerical divergence (when the denominator is within machine epsilon to zero).   
    Notice the order of parameters matches the ordering of the words in the FuncName.
ApplyPhaseFunc accepts optional arguments BitEncoding and PhaseOverrides.
ApplyPhaseFunc[... PhaseOverrides -> rules] first consults whether a basis state's index is included in the list of rules {index -> phase}, and if present, uses the prescribed phase in lieu of evaluating f[index].
    PhaseOverrides which correspond to divergences of named phase functions will be used, in lieu of the divergence parameter.
    For multi-variable functions, each index must be a tuple.
ApplyPhaseFunc[..., BitEncoding -> \"Unsigned\"] interprets each sub-register state as an unsigned binary integer, in {0, ..., 2^numQubits-1}
ApplyPhaseFunc[..., BitEncoding -> \"TwosComplement\"] interprets each sub-register state as a two's complement signed integer in {-2^(N-1), ..., +2^(N-1)-1}, where N is the number of qubits (including the sign qubit).
See ?BitEncoding and ?PhaseOverrides."
    ApplyPhaseFunc::error = "`1`"
    
    CalcPauliStringMatrix::usage = "CalcPauliStringMatrix[pauliString] returns the numerical matrix of the given real-weighted sum of Pauli tensors. The number of qubits is assumed to be the largest Pauli target. This accepts only sums of Pauli products with unique qubits and floating-point coefficients, and is computed numerically."
    CalcPauliStringMatrix::error = "`1`"
    
    CalcPauliExpressionMatrix::usage = "CalcPauliExpressionMatrix[expr] returns the sparse, analytic matrix given by the symbolic expression of Pauli operators, X, Y, Z, Id. The number of qubits is assumed to be the largest Pauli target. Accepts the same inputs as SimplfyPaulis[], and is computed symbolically.
CalcPauliExpressionMatrix[expr, numQb] overrides the assumed number of qubits."
    CalcPauliExpressionMatrix::error = "`1`"
    
    CalcPauliStringMinEigVal::usage = "CalcPauliStringMinEigVal[pauliString] returns the ground-state energy of the given real-weighted sum of Pauli tensors.
CalcPauliStringMinEigVal[pauliString, MaxIterations -> n] specifies to use at most n iterations in the invoked Arnaldi/Lanczos's method"
    CalcPauliStringMinEigVal::error = "`1`"

    DestroyQureg::usage = "DestroyQureg[qureg] destroys the qureg associated with the given ID. If qureg is a Symbol, it will additionally be cleared."
    DestroyQureg::error = "`1`"
    
    GetAmp::usage = "GetAmp[qureg, index] returns the complex amplitude of the state-vector qureg at the given index, indexing from 0.
GetAmp[qureg, row, col] returns the complex amplitude of the density-matrix qureg at index [row, col], indexing from [0,0]."
    GetAmp::error = "`1`"
    
    SetAmp::usage = "SetAmp[qureg, index, amp] modifies the indexed amplitude of the state-vector qureg to complex number amp.
SetAmp[qureg, row, col, amp] modifies the indexed (row, col) amplitude of the density-matrix qureg to complex number amp"
    SetAmp::error = "`1`"
    
    GetQuregMatrix::usage = "GetQuregMatrix[qureg] returns the state-vector or density matrix associated with the given qureg."
    GetQuregMatrix::error = "`1`"
    
    SetQuregMatrix::usage = "SetQuregMatrix[qureg, matr] modifies qureg, overwriting its statevector or density matrix with that passed."
    SetQuregMatrix::error = "`1`"
    
    SetQuregToPauliString::usage = "SetQuregToPauliString[qureg, pauliString] overwrites the given density matrix to become a dense matrix representation of the given pauli string.
The state is likely no longer a valid density matrix but is useful as a persistent Z-basis representation of the pauli string, to be used in functions like CalcDensityInnerProduct[] and CalcExpecPauliStringDerivs[]."
    SetQuregToPauliString::error = "`1`"
    
    GetPauliStringFromCoeffs::usage = "GetPauliStringFromCoeffs[addr] opens or downloads the file at addr (a string, of a file location or URL), and interprets it as a list of coefficients and Pauli codes, converting this to a symbolic weighted sum of Pauli tensors. Each line of the file is a separate term (a Pauli product), with format {coeff code1 code2 ... codeN} (exclude braces) where the codes are in {0,1,2,3} (indicating a I, X, Y, Z term in the product respectively), for an N-qubit operator. Each line must have N+1 terms (including the real decimal coefficient at the beginning)."
    GetPauliStringFromCoeffs::error = "`1`"
    
    GetRandomPauliString::usage = "GetRandomPauliString[numQubits, numTerms, {minCoeff, maxCoeff}] generates a random Pauli string with unique Pauli tensors.
GetRandomPauliString[numQubits, All, {minCoeff, maxCoeff}] will generate all 4^numQubits unique Pauli tensors.
GetRandomPauliString[numQubits, {minCoeff, maxCoeff}] will generate 4 numQubits^4 unique terms / Pauli tensors, unless this exceeds the maximum of 4^numQubits.
GetRandomPauliString[numQubits] will generate random coefficients in [-1, 1].
All combinations of optional arguments are possible."
    GetRandomPauliString::error = "`1`"
    
    CreateRemoteQuESTEnv::usage = "CreateRemoteQuESTEnv[ip, port1, port2] connects to a remote QuESTlink server at ip, at the given ports, and defines several QuEST functions, returning a link object. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link]."
    CreateRemoteQuESTEnv::error = "`1`"
    
    CreateLocalQuESTEnv::usage = "CreateLocalQuESTEnv[fn] connects to a local 'quest_link' executable, located at fn, running single-CPU QuEST. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link].
CreateLocalQuESTEnv[] connects to a 'quest_link' executable in the working directory."
    CreateLocalQuESTEnv::error = "`1`"
    
    CreateDownloadedQuESTEnv::usage = "CreateDownloadedQuESTEnv[] downloads a precompiled single-CPU QuESTlink binary (specific to your operating system) directly from Github, then locally connects to it. This should be called once, before using the QuESTlink API.
CreateDownloadedQuESTEnv[os] forces downloaded of the binary for operating system 'os', which must one of {Windows, Linux, Unix, MacOS, MacOSX}."
    CreateDownloadedQuESTEnv::error = "`1`"
    
    DestroyQuESTEnv::usage = "DestroyQuESTEnv[link] disconnects from the QuEST link, which may be the remote Igor server or a loca instance, clearing some QuEST function definitions (but not those provided by the QuEST package)."
    DestroyQuESTEnv::error = "`1`"

    SetWeightedQureg::usage = "SetWeightedQureg[fac1, q1, fac2, q2, facOut, qOut] modifies qureg qOut to be (facOut qOut + fac1 q1 + fac2 q2). qOut can be one of q1 an q2, and all factors can be complex.
SetWeightedQureg[fac1, q1, fac2, q2, qOut] modifies qureg qOut to be (fac1 q1 + fac2 q2). qOut can be one of q1 an q2.
SetWeightedQureg[fac1, q1, qOut] modifies qureg qOut to be fac1 * q1. qOut can be q1.
SetWeightedQureg[fac, qOut] modifies qureg qOut to be fac qOut; that is, qOut is scaled by factor fac."
    SetWeightedQureg::error = "`1`"
    
    SimplifyPaulis::usage = "SimplifyPaulis[expr] freezes commutation and analytically simplifies the given expression of Pauli operators, and expands it in the Pauli basis. The input expression can include sums, products, powers and non-commuting products of (subscripted) Id, X, Y and Z operators and other Mathematica symbols (including variables defined as Pauli expressions). 
For example, try SimplifyPaulis[ Subscript[Y,0] (a Subscript[X,0] + b Subscript[Z,0] Subscript[X,1])^3 ].
Be careful of performing algebra with Pauli operators outside of SimplifyPaulis[], since Mathematica may erroneously automatically commute them."
    SimplifyPaulis::error = "`1`"

    DrawCircuit::usage = "DrawCircuit[circuit] generates a circuit diagram. The circuit can contain symbolic parameters.
DrawCircuit[circuit, numQubits] generates a circuit diagram with numQubits, which can be more or less than that inferred from the circuit.
DrawCircuit[{circ1, circ2, ...}] draws the total circuit, divided into the given subcircuits. This is the output format of GetCircuitColumns[].
DrawCircuit[{{t1, circ1}, {t2, circ2}, ...}] draws the total circuit, divided into the given subcircuits, labeled by their scheduled times {t1, t2, ...}. This is the output format of GetCircuitSchedule[].
DrawCircuit[{{t1, A1,A2}, {t2, B1,B2}, ...}] draws the total circuit, divided into subcircuits {A1 A2, B1 B2, ...}, labeled by their scheduled times {t1, t2, ...}. This is the output format of InsertCircuitNoise[].
DrawCircuit accepts optional arguments Compactify, DividerStyle, SubcircuitSpacing, SubcircuitLabels, LabelDrawer and any Graphics option. For example, the fonts can be changed with 'BaseStyle -> {FontFamily -> \"Arial\"}'."
    DrawCircuit::error = "`1`"
    
    DrawCircuitTopology::usage = "DrawCircuitTopology[circuit] generates a graph plot of the qubit connectivity implied by the given circuit. The precise nature of the information plotted depends on the following options.
DrawCircuitTopology accepts optional arguments DistinguishBy, ShowLocalGates, ShowRepetitions to modify the presented graph.
DrawCircuitTopology additionally accepts DistinguishedStyles and all options of Graph[], Show[] and LineLegend[] for customising the plot aesthetic."
    DrawCircuitTopology::error = "`1`"

    CalcCircuitMatrix::usage = "CalcCircuitMatrix[circuit] returns an analytic matrix for the given unitary circuit, which may contain symbolic parameters. The number of qubits is inferred from the circuit indices (0 to maximum specified).
CalcCircuitMatrix[circuit] returns an analytic superoperator for the given non-unitary circuit, expressed as a matrix upon twice as many qubits. The result can be multiplied upon a column-flattened density matrix.
CalcCircuitMatrix[circuit, numQubits] forces the number of present qubits.
CalcCircuitMatrix accepts optional argument AsSuperoperator and AssertValidChannels."
    CalcCircuitMatrix::error = "`1`"
    
    GetCircuitGeneralised::usage = "GetCircuitGeneralised[circuit] returns an equivalent circuit composed only of general unitaries (and Matr operators) and Kraus operators of analytic matrices."
    GetCircuitGeneralised::error = "`1`"
    
    GetCircuitSuperoperator::usage = "GetCircuitSuperoperator[circuit] returns the corresponding superoperator circuit upon doubly-many qubits as per the Choi–Jamiolkowski isomorphism. Decoherence channels become Matr[] superoperators.
GetCircuitSuperoperator[circuit, numQubits] forces the circuit to be assumed size numQubits, so that the output superoperator circuit is of size 2*numQubits.
GetCircuitSuperoperator accepts optional argument AssertValidChannels."
    GetCircuitSuperoperator::error = "`1`"
    
    PlotDensityMatrix::usage = "PlotDensityMatrix[matrix] (accepts id or numeric matrix) plots a component (default is magnitude) of the given matrix as a 3D bar plot.
PlotDensityMatrix[matrix1, matrix2] plots both matrix1 and matrix2 simultaneously, and the latter is intended as a \"reference\" state.
PlotDensityMatrix[matrix, vector] converts the state-vector to a density matrix, and plots.
PlotDensityMatrix accepts optional arguments PlotComponent, BarSpacing and all options for Histogram3D. Customising the colour may require overriding the default ColorFunction.
When two matrices are passed, many options (e.g. ChartStyle) can accept a length-2 list."
    PlotDensityMatrix::error = "`1`"
    
    GetCircuitColumns::usage = "GetCircuitColumns[circuit] divides circuit into sub-circuits of gates on unique qubits (i.e. columns), filled from the left. Flatten the result to restore an equivalent but potentially compacted Circuit."
    GetCircuitColumns::error = "`1`"
    
    GetUnsupportedGates::usage = "GetUnsupportedGates[circuit, spec] returns a list of the gates in circuit which either on non-existent qubits or are not present in or satisfy the gate rules in the device specification. The circuit can contain symbolic parameters, though if it cannot be inferred that the parameter satisfies a gate condition, the gate is assumed unsupported.
GetUnsupportedGates[{circ1, circ2, ...}, spec] returns the unsupported gates in each subcircuit, as separate lists.
GetUnsupportedGates[{{t1, circ1}, {t2, circ2}, ...}, spec] ignores the times in the schedule and returns the unsupported gates in each subcircuit, as separate lists."
    GetUnsupportedGates::error = "`1`"
    
    GetCircuitSchedule::usage = "GetCircuitSchedule[circuit, spec] divides circuit into sub-circuits of simultaneously-applied gates (filled from the left), and assigns each a start-time based on the duration of the slowest gate according to the given device specification. The returned structure is {{t1, sub-circuit1}, {t2, sub-circuit2}, ...}, which can be given directly to DrawCircuit[] or ViewCircuitSchedule[].
GetCircuitSchedule[subcircuits, spec] uses the given division (lists of circuits), assumes the gates in each can be performed simultaneously, and performs the same scheduling.
GetCircuitSchedule accepts optional argument ReplaceAliases.
GetCircuitSchedule will take into consideration gates with durations dependent on their scheduled start time."
    GetCircuitSchedule::error = "`1`"
    
    CheckCircuitSchedule::usage = "CheckCircuitSchedule[{{t1, circ1}, {t2, circ2}, ...}, spec] checks whether the given schedule of sub-circuits is compatible with the device specification, can be made compatible, or else if it prescribes overlapping sub-circuit execution (regardless of targeted qubits). Times and gate parameters can be symbolic. All gates in a sub-circuit are assumed applicable simultaneously, even if they target overlapping qubits.
CheckCircuitSchedule returns False if the (possibly symbolic) times cannot possibly be monotonic, nor admit a sufficient duration for any sub-circuit.
CheckCircuitSchedule returns True if the schedule is valid for any assignment of the times and gate parameters.
CheckCircuitSchedule returns a list of symbolic conditions which must be simultaneously satisfied for the schedule to be valid, if it cannot determine so absolutely. These conditions include constraints of both motonicity and duration.
CheckCircuitSchedule will take into consideration gates with durations dependent on their scheduled start time, and circuit variables."
    CheckCircuitSchedule::error = "`1`"
    
    InsertCircuitNoise::usage = "InsertCircuitNoise[circuit, spec] divides the circuit into scheduled subcircuits, then replaces them with rounds of active and passive noise, according to the given device specification. Scheduling is performed by GetCircuitSchedule[]. The output format is {{t1, active, passive}, ...}, which can be given directly to DrawCircuit[], ViewCircuitSchedule[] or ExtractCircuit[].
InsertCircuitNoise[{circ1, circ2, ...}, spec] uses the given list of sub-circuits (output format of GetCircuitColumns[]), assuming each contain gates which can be simultaneously performed.
InsertCircuitNoise[{{t1, circ1}, {t2, circ2}, ...} assumes the given schedule (output format of GetCircuitSchedule[]) of {t1,t2,...} for the rounds of gates and noise. These times can be symbolic.
InsertCircuitNoise accepts optional argument ReplaceAliases.
InsertCircuitNoise can handle gates with time-dependent noise operators and durations."
    InsertCircuitNoise::error = "`1`"
    
    ExtractCircuit::usage = "ExtractCircuit[] returns the prescribed circuit from the outputs of InsertCircuitNoise[], GetCircuitSchedule[] and GetCircuitColumns[]."
    ExtractCircuit::error = "`1`"
    
    ViewCircuitSchedule::usage = "ViewCircuitSchedule[schedule] displays a table form of the given circuit schedule, as output by InsertCircuitNoise[] or GetCircuitSchedule[].
ViewCircuitSchedule accepts all optional arguments of Grid[], for example 'FrameStyle', and 'BaseStyle -> {FontFamily -> \"CMU Serif\"}'."
    ViewCircuitSchedule::error = "`1`"
    
    ViewDeviceSpec::usage = "ViewDeviceSpec[spec] displays all information about the given device specification in table form.
ViewDeviceSpec accepts all optional arguments of Grid[] (to customise all tables), and Column[] (to customise their placement)."
    ViewDeviceSpec::error = "`1`"
    
    CheckDeviceSpec::usage = "CheckDeviceSpec[spec] checks that the given device specification satisfies a set of validity requirements, returning True if so, otherwise reporting a specific error. This is a useful debugging tool when creating a device specification, though a result of True does not gaurantee the spec is valid."
    CheckDeviceSpec::error = "`1`"
    
    GetCircuitInverse::usage = "GetCircuitInverse[circuit] returns a circuit prescribing the inverse unitary operation of the given circuit."
    GetCircuitInverse::error = "`1`"
    
    SimplifyCircuit::usage = "SimplifyCircuit[circuit] returns an equivalent but simplified circuit."
    SimplifyCircuit::error = "`1`"
    
    GetKnownCircuit::usage = "GetKnownCircuit[\"QFT\", qubits]
GetKnownCircuit[\"Trotter\", hamil, order, reps, time]
    (https://arxiv.org/pdf/math-ph/0506007.pdf)
GetKnownCircuit[\"HardwareEfficientAnsatz\", reps, paramSymbol, qubits]
    (https://arxiv.org/pdf/1704.05018.pdf)
GetKnownCircuit[\"TrotterAnsatz\", hamil, order, reps, paramSymbol]
    (https://arxiv.org/pdf/1507.08969.pdf)
GetKnownCircuit[\"LowDepthAnsatz\", reps, paramSymbol, qubits]
    (https://arxiv.org/pdf/1801.01053.pdf)"
    GetKnownCircuit::error = "`1`"
    
    (* BELOW ARE TEMPORARILY MADE PRIVATE DUE TO INSUFFICIENT TESTING BEFORE v012 RELEASE *)
    (*
    GetCircuitsFromChannel::usage = "GetCircuitsFromChannel[channel] returns a list of all pure, analytic circuits which are admitted as possible errors of the input channel (a circuit including decoherence). Coherent noise channels become unitaries weighted by a non-unitary Fac[] operator, while incoherent noise channels become non-trace-preserving Matr[] operators. The sum of the expected values of the (potentially unnormalised) state-vectors output by the returned circuits is equivalent to the expected value of the input channel.
See GetRandomCircuitFromChannel[] to randomly select one of these circuits, weighted by its probability.
See SampleExpecPauliString[] to sample such circuits in order to efficiently approximate the effect of decoherence on an expectation value."
    GetCircuitsFromChannel::error = "`1`"
    
    GetRandomCircuitFromChannel::usage = "GetCircuitsFromChannel[channel] returns a pure, random circuit from the coherent decomposition of the input channel (a circuit including decoherence), weighted by its probability. The average of the expected values of the circuits returned by this function approaches the expected value of the noise channel.
    See SampleExpecPauliString[] to sample such circuits in order to efficiently approximate the effect of decoherence on an expectation value."
    GetRandomCircuitFromChannel::error = "`1`"
    
    SampleExpecPauliString::usage = "SampleExpecPauliString[initQureg, channel, pauliString, numSamples] estimates the expected value of pauliString under the given channel (a circuit including decoherence) upon the state-vector initQureg, through Monte Carlo sampling. This avoids the quadratically greater memory costs of density-matrix simulation, but may need many samples to be accurate.
SampleExpecPauliString[initQureg, channel, pauliString, All] deterministically samples each channel decomposition once.
SampleExpecPauliString[initQureg, channel, pauliString, numSamples, {workQureg1, workQureg2}] uses the given persistent working registers to avoid their internal creation and destruction.
Use option ShowProgress to monitor the progress of sampling."
    SampleExpecPauliString::error = "`1`"
    *)
    
    SampleClassicalShadow::usage = "SampleClassicalShadow[qureg, numSamples] returns a sequence of pseudorandom measurement bases (X, Y and Z) and their outcomes (as bits) when performed on all qubits of the given input state.
\[Bullet] The output has structure { {bases, outcomes}, ...} where bases is a list of Pauli bases (encoded as 1=X, 2=Y, 3=Z) specified per-qubit, and outcomes are the corresponding classical qubit outcomes (0 or 1).
\[Bullet] Both lists are ordered with least significant qubit (index 0) first.
\[Bullet] The output shadow is useful for efficient experimental estimation of quantum state properties, as per Nat. Phys. 16, 1050–1057 (2020)."
    SampleClassicalShadow::error = "`1`"
    
    CalcExpecPauliProdsFromClassicalShadow::usage = "CalcExpecPauliProdsFromClassicalShadow[shadow, prods] returns a list of expected values of each Pauli product, as prescribed by the given classical shadow (e.g. output from SampleClassicalShadow[]).
CalcExpecPauliProdsFromClassicalShadow[shadow, prods, numBatches] divides the shadow into batches, computes the expected values of each, then returns their medians. This may suppress measurement errors. The default numBatches is 10. 
This is the procedure outlined in Nat. Phys. 16, 1050–1057 (2020)."
    CalcExpecPauliProdsFromClassicalShadow::error = "`1`"
    
    
    (*
     * optional arguments to public functions
     *)
     
    BeginPackage["`Option`"]

    WithBackup::usage = "Optional argument to ApplyCircuit, indicating whether to create a backup during circuit evaluation to restore the input state in case of a circuit error. This incurs additional memory (default True). If the circuit contains no error, this option has no effect besides wasting memory."
    
    ShowProgress::usage = "Optional argument to ApplyCircuit and SampleExpecPauliString, indicating whether to show a progress bar during circuit evaluation (default False). This slows evaluation slightly."
    
    PlotComponent::Usage = "Optional argument to PlotDensityMatrix, to plot the \"Real\", \"Imaginary\" component of the matrix, or its \"Magnitude\" (default)."
    
    Compactify::usage = "Optional argument to DrawCircuit, to specify (True or False) whether to attempt to compactify the circuit (or each subcircuit) by left-filling columns of gates on unique qubits (the result of GetCircuitColumns[]). No compactifying may yield better results for circuits with multi-target gates (which invoke swaps)."
    
    DividerStyle::usage = "Optional argument to DrawCircuit, to style the vertical lines separating subcircuits. Use DividerStyle -> None to draw without dividers, and DividerStyle -> Directive[...] to specify multiple styles properties."
    
    SubcircuitSpacing::usage = "Optional argument to DrawCircuit, to specify the horizontal space inserted between subcircuits."
    
    SubcircuitLabels::usage = "Optional argument to DrawCircuit, specifying the list of labels to display between subcircuits. Use 'None' to skip a label while still drawing the divider (except for the first and last divider). Customise these labels with LabelDrawer."
    
    LabelDrawer::usage = "Optional argument to DrawCircuit, to specify a two-argument function for drawing subcircuit labels. For example, Function[{msg,x},Text[msg,{x,-.5}]]. Use LabelDrawer -> None to show no labels."
    
    ShowLocalGates::usage = "Optional argument to DrawCircuitTopology, to specify (True or False) whether single-qubit gates should be included in the plot (as single-vertex loops)."
    
    ShowRepetitions::usage = "Optional argument to DrawCircuitTopology, to specify (True or False) whether repeated instances of gates (or other groups as set by DistinguishBy) in the circuit should each yield a distinct edge.
For example, if ShowRepetitions -> True and DistinguishBy -> \"Qubits\", then a circuit containing three C[Rz] gates between qubits 0 and 1 will produce a graph with three edges between vertices 0 and 1."
    
    DistinguishBy::usage = "Optional argument to DrawCircuitTopology to specify how gates are aggregated into graph edges and legend labels. The possible values (in order of decreasing specificity) are \"Parameters\", \"Qubits\", \"NumberOfQubits\", \"Gates\", \"None\", and a distinct \"Connectivity\" mode.
DistinguishBy -> \"Parameters\" assigns every unique gate (even distinguishing similar operators with different parameters) its own label.
DistinguishBy -> \"Qubits\" discards gate parameters, but respects target qubits, so will assign similar gates (acting on the same qubits) but with different parameters to the same label.
DistinguishBy -> \"NumberOfQubits\" discards gate qubit indices, but respects the number of qubits in a gate. Hence, for example, similar gates controlled on different pairs of qubits will be merged together, but not with the same gate controlled on three qubits.
DistinguishBy -> \"Gates\" respects only the gate type (and whether it is controlled or not), and discards all qubit and parameter information. Hence similar gates acting on different numbers of qubits will be merged to one label. This does not apply to pauli-gadget gates R, which remain distinguished for unique pauli sequences (though discarding qubit indices).
DistinguishBy -> \"None\" performs no labelling or distinguishing of edges.
DistinguishBy -> \"Connectivity\" merges all gates, regardless of type, acting upon the same set of qubits (orderless)."

    DistinguishedStyles::usage = "Optional argument to DrawCircuitTopology, to specify the colours/styles used for each distinguished group (hence ultimately, the edge and legend styles). This must be a list of graphic directives, and will be repeated if it contains too few elements.
DistinguishedStyles -> Automatic will colour the groups by sampling ColorData[\"Rainbow\"]."
        
    ReplaceAliases::usage = "Optional argument to GetCircuitSchedule and InsertCircuitNoise, specifying (True or False) whether to substitute the device specification's alias operators in the output (including in gates and active/passive noise). 
This is False by default, but must be True to pass the output circuits to (for example) ApplyCircuit which don't recognise the alias.
Note if ReplaceAliases -> True, then the output of GetCircuitSchedule might not be compatible as an input to InsertCircuitNoise."

    PhaseOverrides::usage = "Optional argument to ApplyPhaseFunc, specifying overriding values of the phase function for specific sub-register basis states. This can be used to avoid divergences in the phase function.
For a single-variable phase function, like ApplyPhaseFunc[..., x^2, x], the option must have form PhaseOverrides -> {integer -> real, ...}, where 'integer' is the basis state index of which to override the phase.
For a multi-variable phase function, like ApplyPhaseFunc[..., x^2+y^2, {x,y}], the option must have form PhaseOverrides -> {{integer, integer} -> real, ...}, where the basis state indices are specified as a tuple {x,y}.
Note under BitEncoding -> \"TwosComplement\", basis state indices can be negative."
    
    BitEncoding::usage = "Optional argument to ApplyPhaseFunc, specifying how the values of sub-register basis states are encoded in (qu)bits.
BitEncoding -> \"Unsigned\" (default) interprets basis states as natural numbers {0, ..., 2^numQubits-1}.
BitEncoding -> \"TwosComplement\" interprets basis states as two's complement signed numbers, {0, ... 2^(numQubits-1)-1} and {-1, -2, ... -2^(numQubits-1)}. The last qubit in a sub-register list is assumed the sign bit."

    AsSuperoperator::usage = "Optional argument to CalcCircuitMatrix (default Automatic), specifying whether the output should be a 2^N by 2^N unitary matrix (False), or a 2^2N by 2^2N superoperator matrix (True). The latter can capture decoherence, and be multiplied upon column-flattened 2^2N vectors."
    
    AssertValidChannels::usage = "Optional argument to CalcCircuitMatrix and GetCircuitSuperoperator (default True), specifying whether to simplify their outputs by asserting that all channels therein are completely-positive and trace-preserving. For example, this asserts that the argument to a damping channel lies between 0 and 1."
    
    EndPackage[]
    
    
    
    (* 
     * gate symbols    
     *)
     
    BeginPackage["`Gate`"]

    H::usage = "H is the Hadamard gate."
    Protect[H]
        
    X::usage = "X is the Pauli X gate, a.k.a NOT or bit-flip gate."
    Protect[X]
    
    Y::usage = "Y is the Pauli Y gate."
    Protect[Y]
    
    Z::usage = "Z is the Pauli Z gate."
    Protect[Z]
    
    Rx::usage = "Rx[\[Theta]] is a rotation of \[Theta] around the x-axis of the Bloch sphere, Exp[-\[ImaginaryI] \[Theta]/2 X \[CircleTimes] X \[CircleTimes]...]."        
    Protect[Rx]
    
    Ry::usage = "Ry[\[Theta]] is a rotation of \[Theta] around the y-axis of the Bloch sphere, Exp[-\[ImaginaryI] \[Theta]/2 Y \[CircleTimes] Y \[CircleTimes]...]." 
    Protect[Ry]
    
    Rz::usage = "Rz[\[Theta]] is a rotation of \[Theta] around the z-axis of the Bloch sphere, Exp[-\[ImaginaryI] \[Theta]/2 Z \[CircleTimes] Z \[CircleTimes]...]." 
    Protect[Rz]
    
    R::usage = "R[\[Theta], paulis] is the unitary Exp[-\[ImaginaryI] \[Theta]/2 \[CircleTimes] paulis]."   
    Protect[R]
    
    S::usage = "S is the S gate, a.k.a. PI/2 gate."
    Protect[S]
    
    T::usage = "T is the T gate, a.k.a PI/4 gate."
    Protect[T]
    
    U::usage = "U[matrix] is a general unitary gate with any number of target qubits, specified as a unitary square complex matrix. If the matrix breaks normalisation but is intended unitarity, use UNonNorm. To specify a general non-unitary matrix, use Matr."
    Protect[U]
    
    Deph::usage = "Deph[prob] is a 1 or 2 qubit dephasing with probability prob of error."
    Protect[Deph]
    
    Depol::usage = "Depol[prob] is a 1 or 2 qubit depolarising with probability prob of error."
    Protect[Depol]
    
    Damp::usage = "Damp[prob] is 1 qubit amplitude damping with the given decay probability."
    Protect[Damp]
    
    SWAP::usage = "SWAP is a 2 qubit gate which swaps the state of two qubits."
    Protect[SWAP]
    
    M::usage = "M is a destructive measurement gate which measures the indicated qubits in the Z basis. Targeting multiple qubits is the same as applying M to each in-turn, though their outcomes will be grouped in the output of ApplyCircit[]."
    Protect[M]
    
    P::usage = "P[val] is a (normalised) projector onto {0,1} (i.e. a forced measurement) such that the target qubits represent integer val in binary (right most target takes the least significant digit in val).
P[outcome1, outcome2, ...] is a (normalised) projector onto the given {0,1} outcomes. The left most qubit is set to the left most outcome.
The probability of the forced measurement outcome (if hypothetically not forced) is included in the output of ApplyCircuit[]."
    Protect[P]
    
    Kraus::usage = "Kraus[ops] applies a one or two-qubit Kraus map (given as a list of Kraus operators) to a density matrix."
    Protect[Kraus]
    
    KrausNonTP::usage = "KrausNonTP[ops] is equivalent to Kraus[ops] but does not explicitly check that the map is trace-presering. It is still assumed a completely-positive trace-preserving map for internal algorithms, but will tolerate numerical imperfection."
    Protect[KrausNonTP]
    
    G::usage = "G[\[Theta]] applies a global phase rotation of phi, by premultiplying Exp[\[ImaginaryI] \[Theta]]."
    Protect[G]
    
    Id::usage = "Id is an identity gate which effects no change, but can be used for forcing gate alignment in DrawCircuit, or as an alternative to removing gates in ApplyCircuit."
    Protect[Id]
 
    Ph::usage = "Ph is the phase shift gate, which introduces phase factor exp(i*theta) upon state |1...1> of the target and control qubits. The gate is the same under different orderings of qubits, and division between control and target qubits."
    Protect[Ph]
    
    UNonNorm::usage = "UNonNorm[matr] is treated like a general unitary gate U, but with relaxed normalisation conditions on the matrix. This is distinct to gate Matr, which will be internally assumed non-unitary."
    Protect[UNonNorm]
    
    Matr::usage = "Matr[matrix] is an arbitrary operator with any number of target qubits, specified as a completely general (even non-unitary) square complex matrix. Unlike UNonNorm, the given matrix is not internally assumed unitary. It is hence only left-multiplied onto density matrices."
    Protect[Matr]
    
    Fac::usage = "Fac[scalar] is a non-physical operator which multiplies the given complex scalar onto every amplitude of the quantum state. This is directly multiplied onto state-vectors and density-matrices, and may break state normalisation."
    
    (* overriding Mathematica's doc for C[i] as i-th default constant *)
    C::usage = "C is a declaration of control qubits (subscript), which can wrap other gates to conditionally/controlled apply them."
    Protect[C]
 
    EndPackage[]
 
 
 
    (* 
     * device specification keys
     *)
 
    BeginPackage["`DeviceSpec`"]
    
    DeviceDescription::usage = "A description of the specified device."
    
    NumAccessibleQubits::usage = "The number of qubits which can be targeted by user-given circuits in the represented device. These are assumed to be at adjacent indices, starting at 0."
    
    NumTotalQubits::usage = "The number of qubits targeted by all noise in the represented device. This can exceed 'NumAccessibleQubits', since it includes hidden qubits used for advanced noise modelling. Hidden qubits are assumed to start at index 'NumAccessibleQubits'."
    
    Aliases::usage = "Custom aliases for general unitary gates or sub-circuits, recognised by the device specification as elementary gates (optional)."
    
    Gates::usage = "The gates supported by the device, along with their duration and effective operation under active noise."
    
    Qubits::usage = "The qubit properties of the device, such as their passive noise."
    
    NoisyForm::usage = "The channel (expressed as a sub-circuit) describing the noisy, imperfect operation of a device gate."
    
    PassiveNoise::usage = "The channel (expressed as a sub-circuit) describing the passive decoherence of a qubit when not being operated upon by gates."
    
    GateDuration::usage = "The duration of performing a gate on the represented device."
    
    TimeSymbol::usage = "The symbol representing the start time (in a scheduled circuit) of gates and noise channels, which can inform their properties (optional)."
    
    DurationSymbol::usage = "The symbol representing the duration (in a scheduled circuit) of gates or noise channels, which can inform their properties (optional)."
    
    InitVariables::usage = "The function to call at the start of circuit/schedule processing, to re-initialise circuit variables (optional)."
    
    UpdateVariables::usage = "The function to call after each active gate or processed passive noise, to update circuit variables (optional)."
    
    EndPackage[]
    
    
    
    (*
     * deprecated but backwards-compatible API 
     *)
     
    BeginPackage["`Deprecated`"]
    
    CalcExpecPauliProd::usage = "This function is deprecated. Please instead use CalcExpecPauliString."
    CalcExpecPauliSum::usage = "This function is deprecated. Please instead use CalcExpecPauliString."
    ApplyPauliSum::usage = "This function is deprecated. Please instead use ApplyPauliString."
    CalcPauliSumMatrix::usage = "This function is deprecated. Please instead use CalcPauliStringMatrix."
    GetPauliSumFromCoeffs::usage = "This function is deprecated. Please instead use GetPauliStringFromCoeffs."
    MixDamping::usage = "This function is deprecated. Please instead use ApplyCircuit with gate Damp."
    MixDephasing::usage = "This function is deprecated. Please instead use ApplyCircuit with gate Deph."
    MixDepolarising::usage = "This function is deprecated. Please instead use ApplyCircuit with gate Depol."
    MixTwoQubitDephasing::usage = "This function is deprecated. Please instead use ApplyCircuit with gate Deph."
    MixTwoQubitDepolarising::usage = "This function is deprecated. Please instead use ApplyCircuit with gate Depol."
    CalcQuregDerivs::usage = "This function is deprecated. Please instead use ApplyCircuitDerivs."
    
    EndPackage[]
 
 
 
    (* 
     * internal private functions, and definitions of public API
     *)
 
    Begin["`Private`"]
    
    
        
        (*
         * deprecated definitions
         *)
         
        CalcExpecPauliProd[args___] := (
            Message[CalcExpecPauliString::error, "The function CalcExpecPauliProd[] is deprecated. Use CalcExpecPauliString[] or temporarily hide this message using Quiet[]."]; 
            CalcExpecPauliString[args])
        CalcExpecPauliSum[args___] := (
            Message[CalcExpecPauliString::error, "The function CalcExpecPauliSum[] is deprecated. Use CalcExpecPauliString[] or temporarily hide this message using Quiet[]."]; 
            CalcExpecPauliString[args])
        ApplyPauliSum[args___] := (
            Message[ApplyPauliString::error, "The function ApplyPauliSum[] is deprecated, though has still been performed. In future, please use ApplyPauliString[] or temporarily hide this message using Quiet[]."]; 
            ApplyPauliString[args])
        CalcPauliSumMatrix[args___] := (
            Message[CalcPauliStringMatrix::error, "The function CalcPauliSumMatrix[] is deprecated. Use CalcPauliStringMatrix[] or temporarily hide this message using Quiet[]."]; 
            CalcPauliStringMatrix[args])
        GetPauliSumFromCoeffs[args___] := (
            Message[GetPauliStringFromCoeffs::error, "The function GetPauliSumFromCoeffs[] is deprecated. Use GetPauliStringFromCoeffs[] or temporarily hide this message using Quiet[]."]; 
            GetPauliStringFromCoeffs[args])
            
        MixDamping[qureg_Integer, qb_Integer, prob_Real] := (
            Message[ApplyCircuit::error, "The function MixDamping[] is deprecated, though has still been performed. In future, please use ApplyCircuit[] with the Damp[] gate instead, or temporarily hide this message using Quiet[]."];
            ApplyCircuit[qureg, Subscript[Damp,qb][prob]];
            qureg)
        MixDephasing[qureg_Integer, qb_Integer, prob_Real] := (
            Message[ApplyCircuit::error, "The function MixDephasing[] is deprecated, though has still been performed. In future, please use ApplyCircuit[] with the Deph[] gate instead, or temporarily hide this message using Quiet[]."];
            ApplyCircuit[qureg, Subscript[Deph,qb][prob]];
            qureg)
        MixDepolarising[qureg_Integer, qb_Integer, prob_Real] := (
            Message[ApplyCircuit::error, "The function MixDepolarising[] is deprecated, though has still been performed. In future, please use ApplyCircuit[] with the Depol[] gate instead, or temporarily hide this message using Quiet[]."];
            ApplyCircuit[qureg, Subscript[Depol,qb][prob]];
            qureg)
        MixTwoQubitDephasing[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real] := (
            Message[ApplyCircuit::error, "The function MixTwoQubitDephasing[] is deprecated, though has still been performed. In future, please use ApplyCircuit[] with the Deph[] gate instead, or temporarily hide this message using Quiet[]."];
            ApplyCircuit[qureg, Subscript[Deph,qb1,qb2][prob]];
            qureg)
        MixTwoQubitDepolarising[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real] := (
            Message[ApplyCircuit::error, "The function MixTwoQubitDepolarising[] is deprecated, though has still been performed. In future, please use ApplyCircuit[] with the Depol[] gate instead, or temporarily hide this message using Quiet[]."];
            ApplyCircuit[qureg, Subscript[Depol,qb1,qb2][prob]];
            qureg)
            
        CalcQuregDerivs[circuit_, initQureg_, varVals_, derivQuregs_, workQuregs:_:-1] := (
            Message[ApplyCircuitDerivs::error, "The function CalcQuregDerivs[] is deprecated, though has still been attemptedly performed. In future, please use ApplyCircuitDerivs[], or temporarily hide this message using Quiet[]."];
            ApplyCircuitDerivs[initQureg, circuit, varVals, derivQuregs, workQuregs])
            
            
            
        
        (*
         * global convenience functions
         *)
    
        (* report a generic error that the function was passed with bad args (did not evaluate) *)
        invalidArgError[func_Symbol] := (
            Message[func::error, "Invalid arguments. See ?" <> ToString[func]];
            $Failed)
               
        (* opcodes which correlate with the global IDs in circuits.hpp *)
        getOpCode[gate_] :=
	        gate /. {H->0,X->1,Y->2,Z->3,Rx->4,Ry->5,Rz->6,R->7,S->8,T->9,U->10,Deph->11,Depol->12,Damp->13,SWAP->14,M->15,P->16,Kraus->17,G->18,Id->19,Ph->20,KrausNonTP->21,Matr->22,UNonNorm->23,Fac->24,_->-1}
        
        
        
        (*
         * encoding Pauli strings
         *)
         
        pauliCodePatt = X|Y|Z|Id;
        pauliOpPatt = Subscript[pauliCodePatt, _Integer];
        pauliTensorPatt = pauliOpPatt | Verbatim[Times][ Repeated[_?Internal`RealValuedNumericQ,{0,1}], pauliOpPatt.. ];
        
        areUniqueQubits[qubits_List] :=
            CountDistinct[qubits] === Length[qubits]
        
        isValidPauliString[expr_] := Switch[expr,
            pauliOpPatt, 
                True,
            pauliTensorPatt,
                areUniqueQubits[Cases[expr, Subscript[_,q_Integer]:>q]],
            _Plus,
                AllTrue[expr, (MatchQ[#, pauliTensorPatt]&) ]]

        (* X1 *)
        getEncodedPauliString[ Subscript[op:pauliCodePatt, q_Integer] ] := 
            {{1}, {getOpCode@op}, {q}, {1}}
        (* .1 X1 *)
        getEncodedPauliString[ Verbatim[Times][c:_?NumericQ, p:pauliOpPatt.. ] ] := 
            {{c}, getOpCode /@ {p}[[All,1]], {p}[[All,2]], {Length@{p}[[All,2]]}}
        (* X1 X2 *)
        getEncodedPauliString[ Verbatim[Times][p:pauliOpPatt.. ] ] :=
            {{1}, getOpCode /@ {p}[[All,1]], {p}[[All,2]], {Length@{p}[[All,2]]}}
        (* .1 X1 X2 *)
        getEncodedPauliString[ p:pauliTensorPatt ] :=
            {p[[1]], getOpCode /@ Rest[List@@p][[All,1]], Rest[List@@p][[All,2]], Length[p]-1}
        (* .5 X1 X2 + X1 X2 + X1 + .5 X1 *)
        getEncodedPauliString[ s_Plus ] /; AllTrue[List@@s, MatchQ[pauliTensorPatt]] :=
            Join @@@ Transpose[getEncodedPauliString /@ (List @@ s)]
        (* 0.` X1 ... *)
        getEncodedPauliString[ s:Verbatim[Plus][ 0.`, pauliTensorPatt..] ] :=
            getEncodedPauliString @ s[[2;;]]
        
        
        
        (*
         * encoding circuits
         *)
         
        (* checking a product is a valid operator *)
        SetAttributes[isOperatorFormat, HoldAll]
        isOperatorFormat[op_Times] := isCircuitFormat[ReleaseHold[List @@@ Hold[op]]]
        isOperatorFormat[___] := False
        
        (* convert an operator into a circuit spec without commuting gates *)
        SetAttributes[Circuit, HoldAll]
        SetAttributes[Operator, HoldAll]
        Circuit[gate_?isGateFormat] := 
            {gate}
        Circuit[op_?isOperatorFormat] := 
            ReleaseHold[List @@@ Hold[op]]
        Operator[op_?isGateFormat] :=
            {gate}
        Operator[op_?isOperatorFormat] :=
            Reverse @ Circuit @ op
        
        (* convert MMA matrix to a flat format which can be embedded in the circuit param list *)
        codifyMatrix[matr_] :=
            Riffle[Re @ N @ Flatten @ matr, Im @ N @ Flatten @ matr]
            
        (* convert multiple MMA matrices into {#matrices, ... flattened matrices ...} *)
        codifyMatrices[matrs_] :=
            Prepend[Join @@ (codifyMatrix /@ matrs), Length @ matrs]
        
        (* recognising and codifying gates into {opcode, ctrls, targs, params} *)
        gatePatterns = {
            Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[g:U|Matr|UNonNorm,  (targs:__Integer)|{targs:__Integer}][matr:_List]] :> 
                {getOpCode[g], {ctrls}, {targs}, codifyMatrix[matr]},
        	Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]] :> 
                {getOpCode[gate], {ctrls}, {targs}, {args}},
        	Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]] :> 
                {getOpCode[gate], {ctrls}, {targs}, {}},
            Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][R[param_, ({paulis:pauliOpPatt..}|Verbatim[Times][paulis:pauliOpPatt..]|paulis:pauliOpPatt__)]] :>
                {getOpCode[R], {ctrls}, {paulis}[[All,2]], Join[{param}, getOpCode /@ {paulis}[[All,1]]]},
            R[param_, ({paulis:pauliOpPatt..}|Verbatim[Times][paulis:pauliOpPatt..]|paulis:pauliOpPatt__)] :>
                {getOpCode[R], {}, {paulis}[[All,2]], Join[{param}, getOpCode /@ {paulis}[[All,1]]]},
        	Subscript[g:U|Matr|UNonNorm, (targs:__Integer)|{targs:__Integer}][matr:_List] :> 
                {getOpCode[g], {}, {targs}, codifyMatrix[matr]},
            Subscript[Kraus, (targs:__Integer)|{targs:__Integer}][matrs_List] :>
                {getOpCode[Kraus], {}, {targs}, codifyMatrices[matrs]},
            Subscript[KrausNonTP, (targs:__Integer)|{targs:__Integer}][matrs_List] :>
                {getOpCode[KrausNonTP], {}, {targs}, codifyMatrices[matrs]},
            Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__] :> 
                {getOpCode[gate], {}, {targs}, {args}},
        	Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}] :> 
                {getOpCode[gate], {}, {targs}, {}},
            G[arg_] :> 
                {getOpCode[G], {}, {}, {arg}},
            Fac[arg_] :>
                {getOpCode[Fac], {}, {}, {Re @ N @ arg, Re @ Im @ arg}}
        };

        (* converting gate sequence to code lists: {opcodes, ctrls, targs, params} *)
        codifyCircuit[circuit_List] :=
        	circuit /. gatePatterns // Transpose
        codifyCircuit[circuit_] :=
            codifyCircuit @ {circuit}
            
        (* checking circuit format *)
        isGateFormat[Subscript[_Symbol, (__Integer)|{__Integer}]] := True
        isGateFormat[Subscript[_Symbol, (__Integer)|{__Integer}][__]] := True
        isGateFormat[R[_, (pauliOpPatt|{pauliOpPatt..}|Verbatim[Times][pauliOpPatt..])]] := True
        isGateFormat[(G|Fac)[_]] := True
        isGateFormat[___] := False
        isCircuitFormat[circ_List] := AllTrue[circ,isGateFormat]
        isCircuitFormat[circ_?isGateFormat] := True
        isCircuitFormat[___] := False
        
        unpackEncodedCircuit[codes_List] :=
            Sequence[
                codes[[1]], 
                Flatten @ codes[[2]], Length /@ codes[[2]], 
                Flatten @ codes[[3]], Length /@ codes[[3]],
                Flatten[N /@ codes[[4]]], Length /@ codes[[4]]
            ]
            
        
        
        (*
         * ApplyCircuit[]
         *)
            
        (* declaring optional args to ApplyCircuit *)
        Options[ApplyCircuit] = {
            WithBackup -> True,
            ShowProgress -> False
        };
        
        (* applying a sequence of symoblic gates to a qureg. ApplyCircuitInternal provided by WSTP *)
        applyCircuitInner[qureg_, withBackup_, showProgress:0, circCodes__] :=
            ApplyCircuitInternal[qureg, withBackup, showProgress, circCodes]
        applyCircuitInner[qureg_, withBackup_, showProgress:1, circCodes__] :=
            Monitor[
                (* local private variable, updated by backend *)
                calcProgressVar = 0;
                ApplyCircuitInternal[qureg, withBackup, showProgress, circCodes],
                ProgressIndicator[calcProgressVar]
            ]
        ApplyCircuit[qureg_Integer, {}, OptionsPattern[ApplyCircuit]] :=
            {}
        ApplyCircuit[qureg_Integer, circuit_?isCircuitFormat, OptionsPattern[ApplyCircuit]] :=
        	With[
        		{codes = codifyCircuit[circuit]},
        		Which[
                    MemberQ[codes[[1]], -1],
                    Message[ApplyCircuit::error, "Circuit contained an unrecognised gate: " <> ToString@StandardForm@
                        circuit[[ Position[codes[[1]], -1][[1,1]] ]]]; $Failed,
        			Not @ AllTrue[codes[[4]], Internal`RealValuedNumericQ, 2],
                    Message[ApplyCircuit::error, "Circuit contains non-numerical or non-real parameters!"]; $Failed,
                    Not @ Or[OptionValue[WithBackup] === True, OptionValue[WithBackup] === False],
                    Message[ApplyCircuit::error, "Option WithBackup must be True or False."]; $Failed,
                    Not @ Or[OptionValue[ShowProgress] === True, OptionValue[ShowProgress] === False],
                    Message[ApplyCircuit::error, "Option ShowProgress must be True or False."]; $Failed,
                    True,
        			applyCircuitInner[
                        qureg, 
                        If[OptionValue[WithBackup]===True,1,0], 
                        If[OptionValue[ShowProgress]===True,1,0],
                        unpackEncodedCircuit[codes]
                    ]
        		]
        	]
        (* apply a circuit to get an output state without changing input state. CloneQureg provided by WSTP *)
        ApplyCircuit[inQureg_Integer, circuit_?isCircuitFormat, outQureg_Integer, opts:OptionsPattern[ApplyCircuit]] :=
        	Block[{},
        		QuEST`CloneQureg[outQureg, inQureg];
        		ApplyCircuit[outQureg, circuit, opts]
        	]
        ApplyCircuit[inQureg_Integer, {}, outQureg_Integer, opts:OptionsPattern[ApplyCircuit]] := (
            CloneQureg[outQureg, inQureg];
            {}
        )
        (* warnings for old syntax *)
        ApplyCircuit[((_?isCircuitFormat) | {}), _Integer, OptionsPattern[ApplyCircuit]] := (
            Message[ApplyCircuit::error, "As of v0.8, the arguments have swapped order for consistency. Please now use ApplyCircuit[qureg, circuit]."]; 
            $Failed)
        ApplyCircuit[((_?isCircuitFormat) | {}), _Integer, _Integer, OptionsPattern[ApplyCircuit]] := (
            Message[ApplyCircuit::error, "As of v0.8, the arguments have changed order for consistency. Please now use ApplyCircuit[inQureg, circuit, outQureg]."]; 
            $Failed)
        (* error for bad args *)
        ApplyCircuit[___] := invalidArgError[ApplyCircuit]
        
        
        
        (*
         * encoding circuit derivatives
         *)
        
        encodeDerivParams[Subscript[Rx|Ry|Rz|Ph|Damp|Deph|Depol, __][f_], x_] := {D[f,x]}
        encodeDerivParams[R[f_,_], x_] := {D[f,x]}
        encodeDerivParams[G[f_], x_] := {D[f,x]}
        encodeDerivParams[Fac[f_], x_] := With[{df=D[f,x]}, {Re@N@df, Im@N@df}]
        encodeDerivParams[Subscript[U|Matr|UNonNorm, __][matr_], x_] := With[
            {dm = D[matr,x]}, Riffle[Re @ Flatten @ dm, Im @ Flatten @ dm]]
        encodeDerivParams[Subscript[Kraus|KrausNonTP, __][matrs_List], x_] := 
            (Riffle[Re @ Flatten @ #, Im @ Flatten @ #]&) /@ Table[D[m,x] , {m,matrs}]
        encodeDerivParams[Subscript[C, __][g_], x_] := encodeDerivParams[g, x]
        
        encodeDerivCirc[circuit_, varVals_] := Module[{gateInds, varInds, order, encodedCirc, derivParams},

            (* locate the gate indices of the diff variables *)
            gateInds = DeleteDuplicates /@ (Position[circuit, _?(MemberQ[#])][[All, 1]]& /@ varVals[[All,1]]);
            
            (* validate all variables were present in the circuit *)
            If[AnyTrue[gateInds, (# === {} &)],
                Throw @ "One or more variables were not present in the circuit!"];
            
            (* map flat gate indices to relative indices of variables *)
            varInds = Flatten @ Table[ConstantArray[i, Length @ gateInds[[i]]], {i, Length@varVals}];
            gateInds = Flatten @ gateInds;
            
            (* sort info so that gateInds is increasing (for backend optimisations() *)
            order = Ordering[gateInds];
            gateInds = gateInds[[order]];
            varInds = varInds[[order]];
            
            (* encode the circuit for the backend *)
            encodedCirc = codifyCircuit[(circuit /. varVals)];
            
            (* validate that all gates were recognised *)
            If[MemberQ[encodedCirc[[1]], -1],
                Throw["Circuit contained an unrecognised gate: " <> ToString@StandardForm@
                    circuit[[ Position[encodedCirc[[1]], -1][[1,1]] ]]]];
            
            (* validate the circuit contains no unspecified variables *)
            If[Not @ AllTrue[encodedCirc[[4]], Internal`RealValuedNumericQ, 2],
                Throw @ "The circuit contained variables which were not assigned real values."];

            (* differentiate gate args, and pack for backend (without yet making numerical) *)
            derivParams = MapThread[encodeDerivParams, 
                {circuit[[gateInds]], varVals[[varInds,1]]}];
                
            (* validate all gates with diff variables have known derivatives *)
            If[MemberQ[derivParams, encodeDerivParams[_,_]],
                Throw["Cannot differentiate operator " <> 
                    ToString @ StandardForm @ First @ Cases[derivParams, encodeDerivParams[g_,_] :> g] <> "."]];
            
            (* convert packed diff gate args to numerical *)
            derivParams = derivParams /. varVals // N;
            
            (* validate all gate derivatives could be numerically evaluated *)
            If[Not @ AllTrue[Flatten @ derivParams, NumericQ],
                Throw @ "The circuit contained gate derivatives with parameters which could not be numerically evaluated."];
            
            (* return *)
            {encodedCirc, {gateInds, varInds, derivParams}}]
            
        unpackEncodedDerivCircTerms[{gateInds_, varInds_, derivParams_}] :=
            Sequence[gateInds-1, varInds-1, Flatten @ derivParams, Length /@ Flatten /@ derivParams]
            
            
            
        (*
         * derivatives
         *)
         
        ApplyCircuitDerivs[inQureg_Integer, circuit_?isCircuitFormat, vars:{(_ -> _?Internal`RealValuedNumericQ) ..}, outQuregs:{__Integer}, workQuregs:(_Integer|{__Integer}):-1] :=
            Module[
                {ret, encodedCirc, encodedDerivTerms},
                (* check each var corresponds to an out qureg *)
                If[Length[vars] =!= Length[outQuregs],
                    Message[ApplyCircuitDerivs::error, "An equal number of variables and ouptut quregs must be passed."]; Return@$Failed];
                (* encode deriv circuit for backend, throwing any parsing errors *)
                ret = Catch @ encodeDerivCirc[circuit, vars];
                If[Head@ret === String,
                    Message[ApplyCircuitDerivs::error, ret]; Return @ $Failed];
                (* dispatch states, circuit and derivative circuit to backlend *)
                {encodedCirc, encodedDerivTerms} = ret;
                ApplyCircuitDerivsInternal[
                    inQureg, First@{Sequence@@workQuregs}, outQuregs, 
                    unpackEncodedCircuit @ encodedCirc, 
                    unpackEncodedDerivCircTerms @ encodedDerivTerms]]
                    
        ApplyCircuitDerivs[___] := invalidArgError[ApplyCircuitDerivs]  
        
        CalcExpecPauliStringDerivs[initQureg_Integer, circuit_?isCircuitFormat, varVals:{(_ -> _?Internal`RealValuedNumericQ) ..}, paulis_?isValidPauliString, workQuregs:{___Integer}:{}] :=
            Module[
                {ret, encodedCirc, encodedDerivTerms},
                (* encode deriv circuit for backend, throwing any parsing errors *)
                ret = Catch @ encodeDerivCirc[circuit, varVals];
                If[Head@ret === String,
                    Message[CalcExpecPauliStringDerivs::error, ret]; Return @ $Failed];
                (* send to backend, mapping Mathematica indices to C++ indices *)
                {encodedCirc, encodedDerivTerms} = ret;
                CalcExpecPauliStringDerivsInternal[
                    initQureg, workQuregs,
                    unpackEncodedCircuit @ encodedCirc, 
                    unpackEncodedDerivCircTerms @ encodedDerivTerms,
                    Sequence @@ getEncodedPauliString[paulis]]]

        CalcExpecPauliStringDerivs[initQureg_Integer, circuit_?isCircuitFormat, varVals:{(_ -> _?Internal`RealValuedNumericQ) ..}, hamilQureg_Integer, workQuregs:{___Integer}:{}] :=
            Module[
                {ret, encodedCirc, encodedDerivTerms},
                (* encode deriv circuit for backend, throwing any parsing errors *)
                ret = Catch @ encodeDerivCirc[circuit, varVals];
                If[Head@ret === String,
                    Message[CalcExpecPauliStringDerivs::error, ret]; Return @ $Failed];
                (* send to backend, mapping Mathematica indices to C++ indices *)
                {encodedCirc, encodedDerivTerms} = ret;
                CalcExpecPauliStringDerivsDenseHamilInternal[
                    initQureg, hamilQureg, workQuregs,
                    unpackEncodedCircuit @ encodedCirc, 
                    unpackEncodedDerivCircTerms @ encodedDerivTerms]]
            
        CalcExpecPauliStringDerivs[___] := invalidArgError[CalcExpecPauliStringDerivs]
        
        CalcMetricTensor[initQureg_Integer, circuit_?isCircuitFormat, varVals:{(_ -> _?Internal`RealValuedNumericQ) ..}, workQuregs:{___Integer}:{}] :=
            Module[
                {ret, encodedCirc, encodedDerivTerms, retArrs},
                (* encode deriv circuit for backend, throwing any parsing errors *)
                ret = Catch @ encodeDerivCirc[circuit, varVals];
                If[Head@ret === String,
                    Message[CalcMetricTensor::error, ret]; Return @ $Failed];
                (* send to backend, mapping Mathematica indices to C++ indices *)
                {encodedCirc, encodedDerivTerms} = ret;
                data = CalcMetricTensorInternal[
                    initQureg, workQuregs,
                    unpackEncodedCircuit @ encodedCirc, 
                    unpackEncodedDerivCircTerms @ encodedDerivTerms];
                (* reformat output to complex matrix *)
                If[data === $Failed, data, ArrayReshape[
                    MapThread[#1 + I #2 &, {data[[1]], data[[2]]}], 
                    Length[varVals] {1,1}]]]
                    
        CalcMetricTensor[__] := invalidArgError[CalcMetricTensor]
        
        
        
        (*
         * inner products 
         *)
            
        (* compute a matrix of inner products; this can be used in tandem with ApplyCircuitDerivs to populate the Li matrix *)
        CalcInnerProducts[quregIds:{__Integer}] := 
            With[
                {data=CalcInnerProductsMatrixInternal[quregIds],
                len=Length[quregIds]},
                ArrayReshape[
                    MapThread[#1 + I #2 &, {data[[1]], data[[2]]}], 
                    {len, len}
                ]
            ]
        (* computes a vector of inner products <braId|ketIds[i]> *)
        CalcInnerProducts[braId_Integer, ketIds:{__Integer}] := 
            With[
                {data=CalcInnerProductsVectorInternal[braId, ketIds]},
                MapThread[#1 + I #2 &, {data[[1]], data[[2]]}] 
            ]
        (* error for bad args *)
        CalcInnerProducts[___] := invalidArgError[CalcInnerProducts]
            
        (* compute a real symmetric matrix of density inner products *)
        CalcDensityInnerProducts[quregIds:{__Integer}] :=
            ArrayReshape[
                MapThread[
                    #1 + I #2 &,
                    CalcDensityInnerProductsMatrixInternal[quregIds]],
                {Length @ quregIds, Length @ quregIds}
            ]
        (* compute a real vector of density innere products *)
        CalcDensityInnerProducts[rhoId_Integer, omegaIds:{__Integer}] :=
            MapThread[
                #1 + I #2 &,
                CalcDensityInnerProductsVectorInternal[rhoId, omegaIds]]
        (* error for bad args *)
        CalcDensityInnerProducts[___] := invalidArgError[CalcDensityInnerProducts]
        
        
        
        (* 
         * Qureg management 
         *)

        (* destroying a qureg, and clearing the local symbol if recognised *)
        SetAttributes[DestroyQureg, HoldAll];
        DestroyQureg[qureg_Integer] :=
        	DestroyQuregInternal[qureg]
        DestroyQureg[qureg_Symbol] :=
        	Block[{}, DestroyQuregInternal[ReleaseHold@qureg]; Clear[qureg]]
        DestroyQureg[qureg_] :=
            DestroyQuregInternal @ ReleaseHold @ qureg
        DestroyQureg[___] := invalidArgError[DestroyQureg]

        (* get a local matrix representation of the qureg. GetQuregMatrixInternal provided by WSTP *)
        GetQuregMatrix[qureg_Integer] :=
        	With[{data = GetQuregMatrixInternal[qureg]},
        		Which[
        			Or[data === $Failed, data === $Aborted],
        			data,
        			data[[2]] === 0,
        			MapThread[#1 + I #2 &, {data[[3]], data[[4]]}],
        			data[[2]] === 1,
        			Transpose @ ArrayReshape[
        				MapThread[#1 + I #2 &, {data[[3]], data[[4]]}], 
        				{2^data[[1]],2^data[[1]]}]
        		]
        	]
        GetQuregMatrix[___] := invalidArgError[GetQuregMatrix]

        (* overwrite the state of a qureg. InitStateFromAmps provided by WSTP *)
        SetQuregMatrix[qureg_Integer, elems_List] :=
        	With[{flatelems = N @ 
        		Which[
        			(* vectors in various forms *)
        			Length @ Dimensions @ elems === 1,
        				elems,
        			First @ Dimensions @ elems === 1,
        				First @ elems,
        			(Dimensions @ elems)[[2]] === 1,
        				First @ Transpose @ elems,
        			(* density matrices *)
        			SquareMatrixQ @ elems,
        				Flatten @ Transpose @ elems
        		]},
        		QuEST`InitStateFromAmps[qureg, Re[flatelems], Im[flatelems]]
        	]
        SetQuregMatrix[___] := invalidArgError[SetQuregMatrix]
        
        SetQuregToPauliString[qureg_Integer, hamil_?isValidPauliString] :=
            SetQuregToPauliStringInternal[qureg, Sequence @@ getEncodedPauliString[hamil]]
        SetQuregToPauliString[___] := invalidArgError[SetQuregToPauliString]
        
        
        (*
         * Pauli strings
         *)

        invalidPauliScalarError[caller_] := (
            Message[caller::error, "The Pauli string contains a scalar. Perhaps you meant to multiply it onto an identity (Id) operator."]; 
            $Failed)
            
        CalcExpecPauliString[qureg_Integer, paulis_?isValidPauliString, workspace_Integer] :=
            CalcExpecPauliStringInternal[qureg, workspace, Sequence @@ getEncodedPauliString[paulis]]
        CalcExpecPauliString[_Integer, Verbatim[Plus][_?NumericQ, ___], _Integer] := 
            invalidPauliScalarError[CalcExpecPauliString]
        CalcExpecPauliString[___] := invalidArgError[CalcExpecPauliString]

        ApplyPauliString[inQureg_Integer, paulis_?isValidPauliString, outQureg_Integer] :=
            ApplyPauliStringInternal[inQureg, outQureg, Sequence @@ getEncodedPauliString[paulis]]
        ApplyPauliString[_Integer, Verbatim[Plus][_?NumericQ, ___], _Integer] := 
            invalidPauliScalarError[ApplyPauliString]
        ApplyPauliString[___] := invalidArgError[ApplyPauliString]
        
        getFullHilbertPauliMatrix[numQ_][Subscript[s_,q_]] := Module[
        	{m=ConstantArray[SparseArray @ IdentityMatrix[2], numQ]},
        	m[[q+1]] = SparseArray @ PauliMatrix[s /. {Id->0, X->1,Y->2,Z->3}];
        	If[Length[m]>1, KroneckerProduct @@ (Reverse @ m), First @ m]]
            
        SetAttributes[CalcPauliExpressionMatrix, HoldAll]
        CalcPauliExpressionMatrix[h_, nQb_] := With[
        	{hFlat = SimplifyPaulis[h]},
        	ReleaseHold[
        		HoldForm[hFlat] /. Verbatim[Times][a___, b:pauliOpPatt, c___] :>
        			RuleCondition @ Times[
        				Sequence @@ Cases[{a,b,c}, Except[pauliOpPatt]],
        				Dot @@ getFullHilbertPauliMatrix[nQb] /@ Cases[{a,b,c}, pauliOpPatt]
        			] /. p:pauliOpPatt :> RuleCondition @ getFullHilbertPauliMatrix[nQb][p]]]
        CalcPauliExpressionMatrix[h_] := With[
            {hFlat = SimplifyPaulis[h]},
            {nQb = Max[1 + Cases[{hFlat}, Subscript[(Id|X|Y|Z), q_]:>q, Infinity]]},
            CalcPauliExpressionMatrix[hFlat, nQb]]
        CalcPauliExpressionMatrix[___] := invalidArgError[CalcPauliExpressionMatrix]
        
        CalcPauliStringMinEigVal[paulis_?isValidPauliString, MaxIterations -> its_Integer] := With[
            {matr = CalcPauliExpressionMatrix[paulis]},
            - First @ Eigenvalues[- matr, 1, Method -> {"Arnoldi", MaxIterations -> its, "Criteria" -> "RealPart"}]]
        CalcPauliStringMinEigVal[paulis_?isValidPauliString] :=
            CalcPauliStringMinEigVal[paulis, MaxIterations -> 10^5]
        CalcPauliStringMinEigVal[___] := invalidArgError[CalcPauliStringMinEigVal]
        
        CalcPauliStringMatrix[paulis_?isValidPauliString] := With[
            {pauliCodes = getEncodedPauliString[paulis]},
            {elems = CalcPauliStringMatrixInternal[1+Max@pauliCodes[[3]], Sequence @@ pauliCodes]},
            If[elems === $Failed, elems, 
                (#[[1]] + I #[[2]])& /@ Partition[elems,2] // Transpose]]
        CalcPauliStringMatrix[Verbatim[Plus][_?NumericQ, ___]] :=
            invalidPauliScalarError[CalcPauliStringMatrix]
        CalcPauliStringMatrix[___] := invalidArgError[CalcPauliStringMatrix]
        
        GetPauliStringFromCoeffs[addr_String] :=
            Plus @@ (#[[1]] If[ 
                    AllTrue[ #[[2;;]], PossibleZeroQ ],
                    Subscript[Id, 0],
                    Times @@ MapThread[
                    (   Subscript[Switch[#2, 0, Id, 1, X, 2, Y, 3, Z], #1 - 1] /. 
                        Subscript[Id, _] ->  Sequence[] & ), 
                        {Range @ Length @ #[[2 ;;]], #[[2 ;;]]}
                    ]
                ] &) /@ ReadList[addr, Number, RecordLists -> True];
        GetPauliStringFromCoeffs[___] := invalidArgError[GetPauliStringFromCoeffs]
        
        GetRandomPauliString[
            numQubits_Integer?Positive, numTerms:(_Integer?Positive|Automatic|All):Automatic, 
            {minCoeff_?Internal`RealValuedNumericQ, maxCoeff_?Internal`RealValuedNumericQ}
        ] := With[
            {numUniqueTensors = 4^numQubits},
            (* give warning if too many terms requested *)
            If[ NumericQ[numTerms] && numTerms > numUniqueTensors,
                Message[GetRandomPauliString::error, "More terms were requested than there are unique Pauli tensors. Hide this warning with Quiet[]."]]; 
            With[
                {strings = Table[
                    (* generate uniformly random coefficients *)
                    RandomReal[{minCoeff,maxCoeff}] * 
                    Times @@ (
                        (* generate uniformly random but unique Pauli tensors *)
                        MapThread[Subscript[#1, #2]&, {
                            IntegerDigits[tensorInd, 4, numQubits] /. {0->Id,1->X,2->Y,3->Z},
                            Range[0,numQubits-1]}
                            ] /. Subscript[Id, _]->Nothing /. {} -> {Subscript[Id, 0]}),
                        {tensorInd, RandomSample[0;;(numUniqueTensors-1), 
                    (* potentially override the number of terms/tensors *)
                    Min[numTerms /. {Automatic -> 4 numQubits^4, All -> numUniqueTensors}, numUniqueTensors]]}]},
                (* append an Id with max target qubits on the ened for user convenience, if not already a max target  *)
                Plus @@ If[
                    FreeQ[ Last @ strings, Subscript[_Symbol, numQubits-1]],
                    Append[Most @ strings, (Last @ strings) Subscript[Id,numQubits-1]],
                    strings]]]
        GetRandomPauliString[numQubits_Integer?Positive, numTerms:(_Integer?Positive|Automatic|All):Automatic] :=
            GetRandomPauliString[numQubits, numTerms, {-1,1}]
        GetRandomPauliString[___] := invalidArgError[GetRandomPauliString]
        
        
        
        (*
         * Analytic and numerical channel decompositions for statevector simulation
         *)
         
        convertOpToPureCircs[Subscript[(Kraus|KrausNonTP), q__][matrs:{ __?MatrixQ }]] := 
        	Circuit /@ Subscript[Matr, q] /@ matrs
        convertOpToPureCircs[Subscript[Deph, q_][x_]] := {
        	Circuit[ Fac@Sqrt[1-x] ],
        	Circuit[ Fac@Sqrt[x] Subscript[Z, q] ] }
        convertOpToPureCircs[Subscript[Deph, q1_,q2_][x_]] := {
        	Circuit[ Fac@Sqrt[1-x] ],
        	Circuit[ Fac@Sqrt[x/3] Subscript[Z, q1] ],
        	Circuit[ Fac@Sqrt[x/3] Subscript[Z, q2] ],
        	Circuit[ Fac@Sqrt[x/3] Subscript[Z, q1] Subscript[Z, q2] ]}
        convertOpToPureCircs[Subscript[Depol, q_][x_]] := {
        	Circuit[ Fac@Sqrt[1-x] ],
        	Circuit[ Fac@Sqrt[x/3] Subscript[X, q] ],
        	Circuit[ Fac@Sqrt[x/3] Subscript[Y, q] ],
        	Circuit[ Fac@Sqrt[x/3] Subscript[Z, q] ]}
        convertOpToPureCircs[Subscript[Depol, q1_,q2_][x_]] := 
        	Join[{{Fac@Sqrt[1-x]}}, Rest @ Flatten[
        		Table[{Fac@Sqrt[x/15], Subscript[a, q1], Subscript[b, q2]}, {a,{Id,X,Y,Z}}, {b,{Id,X,Y,Z}}] /. Subscript[Id, _] -> Nothing, 1]]
        convertOpToPureCircs[g:Subscript[Damp, q_][x_]] :=
        	convertOpToPureCircs @ First @ GetCircuitGeneralised @ g
        convertOpToPureCircs[g_] :=
        	{{g}} (* non-decoherence ops have no alternatives *)
            
        GetCircuitsFromChannel[circ_List /; isCircuitFormat[circ] ] := With[
            {choices = convertOpToPureCircs /@ circ},
            {numCircs = Times @@ Length /@ choices},
            If[Ceiling @ Log2[numCircs] > $SystemWordLength,
                Message[GetCircuitsFromChannel::error, "The number of unique circuit decompositions exceeds 2^$SystemWordLength and cannot be enumerated."]; $Failed,
        	       Flatten /@ Tuples[choices]]]
        GetCircuitsFromChannel[gate_?isCircuitFormat] :=
            GetCircuitsFromChannel @ {gate}
        GetCircuitsFromChannel[___] := invalidArgError[GetCircuitsFromChannel]
        
        GetRandomCircuitFromChannel[ channel_List /; isCircuitFormat[channel] ] := With[

        	(* get circuit decompositions of each  operator *)
        	{circs = convertOpToPureCircs /@ channel},
        	
        	(* infer the probabilities from Fac[] in coherent noise, and assert uniform incoherent noise *)
        	{probs = Table[
        		Times @@ (Cases[choice, Fac[x_]:>Abs[x]^2] /. {}->{1/N@Length@choices}), 
        		{choices, circs}, {choice, choices}]},
        	
        	(* validate probabilities *)
        	If[ Not[And @@ Table[
        		VectorQ[probset, Internal`RealValuedNumericQ] &&
        		AllTrue[probset, (0 <= # <= 1&)] &&
        		Abs[Total[probset] - 1] < 10^6 $MachineEpsilon, 
        		{probset, probs}]],
        			Message[GetRandomCircuitFromChannel::error, "The probabilities of a decomposition of a decoherence operator were invalid and/or unnormalised."];
        			Return[$Failed]];

        	(* randomly select a pure circuit from each channel decomposition *)
        	Flatten @ MapThread[
        		With[{choice=RandomChoice[#1 -> #2]}, 
        			(* we must multiply the matrices of incoherent noise with the asserted uniform probability *)
        			If[ Length[#1]>1 && Not @ MemberQ[choice, Fac[_]],
        				(* which we can equivalently perform with a Fac gate *)
        				Join[{Fac[Sqrt@N@Length[#1]], choice}],
        				choice /. Fac[_]->Nothing]] &,
        		{probs, circs}]
        ]
        GetRandomCircuitFromChannel[operator_?isCircuitFormat] :=
            GetRandomCircuitFromChannel @ {operator}
        GetRandomCircuitFromChannel[___] := invalidArgError[GetRandomCircuitFromChannel]
        
        Options[SampleExpecPauliString] = {
            ShowProgress -> False
        };
        
        sampleExpecPauliStringInner[True, args__] :=
            Monitor[
                (* local private variable, updated by backend *)
                calcProgressVar = 0;
                SampleExpecPauliStringInternal[1, args],
                ProgressIndicator[calcProgressVar]]
        sampleExpecPauliStringInner[False, args__] :=
            SampleExpecPauliStringInternal[0, args]
         
        SampleExpecPauliString[qureg_Integer, channel_?isCircuitFormat, paulis_?isValidPauliString, numSamples:(_Integer|All), {work1_Integer, work2_Integer}, OptionsPattern[]] /; (work1 === work2 === -1 || And[work1 =!= -1, work2 =!= -1]) :=
            If[numSamples =!= All && numSamples >= 2^63, 
                Message[SampleExpecPauliString::error, "The requested number of samples is too large, and exceeds the maximum C long integer (2^63)."]; $Failed,
                With[{codes = codifyCircuit[channel]},
                    If[
                        Not @ AllTrue[codes[[4]], Internal`RealValuedNumericQ, 2],
                        Message[SampleExpecPauliString::error, "Circuit contains non-numerical or non-real parameters!"]; $Failed,
                        sampleExpecPauliStringInner[
                            OptionValue[ShowProgress],
                            qureg, work1, work2, numSamples /. (All -> -1),
                            unpackEncodedCircuit[codes],
                            Sequence @@ getEncodedPauliString[paulis]]]]]
        SampleExpecPauliString[qureg_Integer, channel_?isCircuitFormat, paulis_?isValidPauliString, numSamples:(_Integer|All), opts:OptionsPattern[]] :=
            SampleExpecPauliString[qureg, channel, paulis, numSamples, {-1, -1}, opts]
        SampleExpecPauliString[___] := invalidArgError[SampleExpecPauliString]
        
        SampleClassicalShadow[qureg_Integer, numSamples_Integer] /; (numSamples >= 2^63) := (
            Message[SampleClassicalShadow::error, "The requested number of samples is too large, and exceeds the maximum C long integer (2^63)."];
            $Failed)
        SampleClassicalShadow[qureg_Integer, numSamples_Integer] := 
            With[
                {data = SampleClassicalShadowStateInternal[qureg, numSamples]},
                If[data === $Failed, data, 
                    Transpose[{
                        Partition[ data[[2]], data[[1]]],
                        Partition[ data[[3]], data[[1]]]}]]]
        SampleClassicalShadow[___] := invalidArgError[SampleClassicalShadow]
    
        CalcExpecPauliProdsFromClassicalShadow[shadow_List, prods:{__:pauliTensorPatt}, numBatches_Integer:10] := 
            If[
                Not @ MatchQ[Dimensions[shadow], {nSamps_, 2, nQb_}],
                (Message[CalcExpecPauliProdsFromClassicalShadow::error, "The classical shadow input must be a list " <>
                    "(length equal to the number of samples) of length-2 sublists, each of length equal to the number 
                    of qubits. This is the format {{{bases,outcomes}},...}, matching that output by SampleClassicalShadow[]."];
                    $Failed),
                With[
                    {ops = (List @@@ prods) /. {p_:pauliCodePatt, q_Integer} :> {Subscript[p, q]}},
                    {numQb = Length @ First @ First @ shadow,
                     prodPaulis = ops[[All, All, 1]] /. {X->1,Y->2,Z->3},
                     prodQubits = ops[[All, All, 2]]},
                    CalcExpecPauliProdsFromClassicalShadowInternal[
                        numQb, numBatches, Length[shadow], Flatten @ shadow[[All,1]], Flatten @ shadow[[All,2]],
                        Flatten @ prodPaulis, Flatten @ prodQubits, Length /@ prodPaulis]
                ]
            ]    
         CalcExpecPauliProdsFromClassicalShadow[___] := invalidArgError[CalcExpecPauliProdsFromClassicalShadow]

    
        
        (*
         * QuESTEnv management 
         *)
        
        getIgorLink[id_] :=
        	LinkConnect[
        		With[{host="@129.67.85.74",startport=50000},
        		ToString[startport+id] <> host <> "," <> ToString[startport+id+100] <> host],
        		LinkProtocol->"TCPIP"]
                
        getRemoteLink[ip_, port1_, port2_] :=
        	LinkConnect[
        		With[{host="@"<>ip},
        		ToString[port1] <> host <> "," <> ToString[port2] <> host],
        		LinkProtocol->"TCPIP"]
                    
        CreateRemoteQuESTEnv[ip_String, port1_Integer, port2_Integer] := Install @ getRemoteLink[ip, port1, port2]
        CreateRemoteQuESTEnv[___] := invalidArgError[CreateRemoteQuESTEnv]
                    
        CreateLocalQuESTEnv[arg_:"quest_link"] := With[
            {fn = arg <> If[$OperatingSystem === "Windows", ".exe", ""]},
            If[
                FileExistsQ[fn], 
                Install[fn],  
                Message[CreateLocalQuESTEnv::error, "Local quest_link executable not found!"]; $Failed]
            ]
        CreateLocalQuESTEnv[__] := invalidArgError[CreateLocalQuESTEnv] (* no args is valid *)
            
        getExecFn["MacOS"|"MacOSX"] = "macos_x86_quest_link";
        getExecFn["MacOS M1"|"MacOSX M1"] = "macos_arm_quest_link";
        getExecFn["Windows"] = "windows_quest_link.exe";
        getExecFn["Linux"|"Unix"] = "linux_quest_link";
        CreateDownloadedQuESTEnv[os:("MacOS"|"MacOSX"|"MacOS M1"|"MacOSX M1"|"Windows"|"Linux"|"Unix")] := 
            Module[{url,resp,log,fn,exec},
                (* attempt download *)
                log = PrintTemporary["Downloading..."];
                url = "https://github.com/QTechTheory/QuESTlink/raw/main/Binaries/" <> getExecFn[os];
                fn = FileNameJoin[{Directory[], "quest_link"}];
                resp = URLDownload[url, fn,  {"File", "StatusCode"}, TimeConstraint->10];
                NotebookDelete[log];
                (* check response *)
                If[resp["StatusCode"] >= 400,
                    Message[CreateDownloadedQuESTEnv::error, "Download failed; returned status code " <> ToString @ resp["StatusCode"]];
                    Return[$Failed, Module]
                ];
                (* install *)
                log = PrintTemporary["Installing..."];
                If[os =!= "Windows", 
                    Run["chmod +x " <> fn];
                    (* esoteric permissions problem. Machine gun hot fix, pow pow *)
                    Run["chmod +x " <> (ToString @@ resp["File"])];
                    Run["chmod +x " <> FileNameTake[fn]];
                ];
                exec = Install @ resp["File"];
                NotebookDelete[log];
                exec
            ]
        CreateDownloadedQuESTEnv[] := Module[
            {os = $OperatingSystem},
            If[ (os === "MacOS" || os === "MacOSX") && Not @ StringContainsQ[$System, "x86"],
                os = "MacOS M1"
            ];
            CreateDownloadedQuESTEnv[os]
        ]
        CreateDownloadedQuESTEnv[__] := (
            Message[CreateDownloadedQuESTEnv::error, "Supported operating systems are Windows, Linux, Unix, MacOS, MacOSX, MacOS M1, MacOSX M1"]; 
            $Failed)
                    
        DestroyQuESTEnv[link_] := Uninstall @ link
        DestroyQuESTEnv[___] := invalidArgError[DestroyQuESTEnv]
        
        
        
        (*
         * state getters and setters 
         *)
         
        (* Im[0.] = 0, how annoying *)
        SetWeightedQureg[fac1_?NumericQ, q1_Integer, fac2_?NumericQ, q2_Integer, facOut_?NumericQ, qOut_Integer] :=
            SetWeightedQuregInternal[
                Re @ N @ fac1, N @ Im @ N @ fac1, q1,
                Re @ N @ fac2, N @ Im @ N @ fac2, q2,
                Re @ N @ facOut, N @ Im @ N @ facOut, qOut
            ]
        SetWeightedQureg[fac1_?NumericQ, q1_Integer, fac2_?NumericQ, q2_Integer, qOut_Integer] :=
            SetWeightedQuregInternal[
                Re @ N @ fac1, N @ Im @ N @ fac1, q1,
                Re @ N @ fac2, N @ Im @ N @ fac2, q2,
                0., 0., qOut
            ]
        SetWeightedQureg[fac1_?NumericQ, q1_Integer, qOut_Integer] :=
            SetWeightedQuregInternal[
                Re @ N @ fac1, N @ Im @ N @ fac1, q1,
                0., 0., q1,
                0., 0., qOut
            ]
        SetWeightedQureg[fac_?NumericQ, qOut_Integer] :=
            SetWeightedQuregInternal[
                0., 0., qOut,
                0., 0., qOut,
                Re @ N @ fac, N @ Im @ N @ fac, qOut
            ]
        SetWeightedQureg[___] := invalidArgError[SetWeightedQureg]
        
        GetAmp[qureg_Integer, index_Integer] := GetAmpInternal[qureg, index, -1]
        GetAmp[qureg_Integer, row_Integer, col_Integer] := GetAmpInternal[qureg, row, col]
        GetAmp[___] := invalidArgError[GetAmp]
        
        SetAmp[qureg_Integer, index_Integer, amp_?NumericQ] := SetAmpInternal[qureg, N@Re@N@amp, N@Im@N@amp, index, -1]
        SetAmp[qureg_Integer, row_Integer, col_Integer, amp_?NumericQ] := SetAmpInternal[qureg, N@Re@N@amp, N@Im@N@amp, row, col]
        SetAmp[___] := invalidArgError[SetAmp]
        
        
        
        (*
         * phase functions
         *)
        
        (* for extracting {coeffs}, {exponents} from a 1D exponential-polynomial *)
        extractCoeffExpo[s_Symbol][c_?NumericQ] := {c,0}
        extractCoeffExpo[s_Symbol][s_Symbol] := {1,1}
        extractCoeffExpo[s_Symbol][Verbatim[Times][c_?NumericQ, s_Symbol]] := {c,1}
        extractCoeffExpo[s_Symbol][Verbatim[Power][s_Symbol,e_?NumericQ]] := {1,e}
        extractCoeffExpo[s_Symbol][Verbatim[Times][c_?NumericQ, Verbatim[Power][s_Symbol,e_?NumericQ]]] := {c,e}
        extractCoeffExpo[s_Symbol][badTerm_] := {$Failed, badTerm}
        extractExpPolyTerms[poly_Plus, s_Symbol] :=
        	extractCoeffExpo[s] /@ List @@ poly
        extractExpPolyTerms[term_, s_Symbol] :=
        	{extractCoeffExpo[s] @ term}
            
        (* for extracting {coeffs}, {exponents} from an n-D exponential-polynomial *)
        extractMultiExpPolyTerms[terms_List, symbs:{__Symbol}] := 
        	Module[{coeffs,powers, cp, badterms},
        		coeffs = Association @@ Table[s->{},{s,symbs}];
        		powers = Association @@ Table[s->{},{s,symbs}];
        		badterms = {};
        		(* for each term... *)
        		Do[
        			(* attempting extraction of term via each symbol *)
        			Do[
        				cp = extractCoeffExpo[s][term];
        				If[ First[cp] =!= $Failed,
        					AppendTo[coeffs[s], cp[[1]]];
        					AppendTo[powers[s], cp[[2]]];
        					Break[]],
        				{s, symbs}];
        			(* if no symbol choice admitted a recognised term, record term *)
        			If[ First[cp] === $Failed,
        				AppendTo[badterms, cp]];
        			(* otherwise, proceed through all terms *)
        			, {term, terms}];
        		(* return bad terms if encountered, else term info *)
        		If[badterms === {}, 
        			{Values @ coeffs, Values @ powers},
        			badterms]]
        extractMultiExpPolyTerms[poly_Plus, symbs:{__Symbol}] := 
            extractMultiExpPolyTerms[List @@ poly, symbs]
        extractMultiExpPolyTerms[term_, symbs:{__Symbol}] := 
            extractMultiExpPolyTerms[{term}, symbs]
            
        bitEncodingFlags = {  (* these must match the values of the enum bitEncoding in QuEST.h *)
            "Unsigned" -> 0,
            "TwosComplement" -> 1
        };
        phaseFuncFlags = {    (* these must match the values of the enum phaseFunc in QuEST.h *)
            "Norm" -> 0,
            "ScaledNorm" -> 1,
            "InverseNorm" -> 2,
            "ScaledInverseNorm" -> 3,
            "ScaledInverseShiftedNorm" -> 4,
            
            "Product" -> 5,
            "ScaledProduct" -> 6,
            "InverseProduct" -> 7,
            "ScaledInverseProduct" -> 8,
            
            "Distance" -> 9,
            "ScaledDistance" -> 10,
            "InverseDistance" -> 11,
            "ScaledInverseDistance" -> 12,
            "ScaledInverseShiftedDistance" -> 13
        };
        Options[ApplyPhaseFunc] = {
            BitEncoding -> "Unsigned",
            PhaseOverrides -> {}
        };
        
        (* single-variable exponential polynomial func *)
        ApplyPhaseFunc[qureg_Integer, reg:{__Integer}, func_, symb_Symbol, OptionsPattern[]] := With[
            {terms = extractExpPolyTerms[N @ func, symb]},
            {badterms = Cases[terms, {$Failed, bad_} :> bad]},
            {overs = OptionValue[PhaseOverrides]},
            Which[
                Length[badterms] > 0,
                    (Message[ApplyPhaseFunc::error, "The phase function, which must be an exponential-polynomial, contained an unrecognised term of the form " <> ToString@StandardForm@First@badterms <> "."]; 
                    $Failed),
                Not @ MemberQ[bitEncodingFlags[[All,1]], OptionValue[BitEncoding]],
                    (Message[ApplyPhaseFunc::error, "Invalid option for BitEncoding. Must be one of " <> ToString@bitEncodingFlags[[All, 1]] <> "."]; 
                    $Failed),
                Not @ MatchQ[overs, {(_Integer -> _?Internal`RealValuedNumericQ) ...}],
                    (Message[ApplyPhaseFunc::error, "Invalid one-dimensional PhaseOverrides, which must be of the form {integer -> real, ...}"]; 
                    $Failed),
                True,
                    ApplyPhaseFuncInternal[qureg, reg, OptionValue[BitEncoding] /. bitEncodingFlags, terms[[All,1]], terms[[All,2]], overs[[All,1]], N @ overs[[All,2]]]]]
            
        (* multi-variable exponential polynomial func *)
        ApplyPhaseFunc[qureg_Integer, regs:{{__Integer}..}, func_, symbs:{__Symbol}, OptionsPattern[]] := With[
            {terms = extractMultiExpPolyTerms[N @ func, symbs]},
            {badterms = Cases[terms, {$Failed, bad_} :> bad]},
            {coeffs = First[terms], exponents=Last[terms]},
            {overs = OptionValue[PhaseOverrides]},
            Which[
                Not @ DuplicateFreeQ @ symbs,
                    (Message[ApplyPhaseFunc::error, "The list of phase function symbols must be unique."];
                    $Failed),
                Length[regs] =!= Length[symbs],
                    (Message[ApplyPhaseFunc::error, "Each delimited sub-register of qubits must correspond to a unique symbol in the phase function."];
                    $Failed),
                Length[badterms] > 0,
                    (Message[ApplyPhaseFunc::error, "The phase function, which must be an exponential-polynomial, contained an unrecognised term of the form " <> ToString@StandardForm@First@badterms <> "."]; 
                     $Failed),
                Not @ MemberQ[bitEncodingFlags[[All,1]], OptionValue[BitEncoding]],
                    (Message[ApplyPhaseFunc::error, "Invalid option for BitEncoding. Must be one of " <> ToString@bitEncodingFlags[[All, 1]] <> "."]; 
                    $Failed),
                Not[ (overs === {}) || And[
                        MatchQ[overs, {({__Integer} -> _?Internal`RealValuedNumericQ) ..}],
                        Equal @@ Length /@ overs[[All,1]], 
                        Length[regs] === Length@overs[[1,1]] ] ],
                    (Message[ApplyPhaseFunc::error, "Invalid PhaseOverrides. Each overriden phase index must be specified as an n-tuple, where n is the number of sub-registers and symbols, pointing to a real number. For example, ApplyPhaseFunc[..., {x,y}, PhaseOverrides -> { {0,0} -> PI, ... }]."];
                     $Failed),
                True,
                    ApplyMultiVarPhaseFuncInternal[qureg, Flatten[regs], Length/@regs, OptionValue[BitEncoding] /. bitEncodingFlags, Flatten[coeffs], Flatten[exponents], Length/@coeffs, Flatten[overs[[All,1]]], N @ overs[[All,2]]]]]

        (* parameterised named func (multi-variable) *)
        ApplyPhaseFunc[qureg_Integer, regs:{{__Integer}..}, {func_String, params___?Internal`RealValuedNumericQ}, OptionsPattern[]] := With[
            {overs = OptionValue[PhaseOverrides]},
            Which[
                Not @ MemberQ[phaseFuncFlags[[All,1]], func],
                    (Message[ApplyPhaseFunc::error, "The named phase function must be one of " <> ToString[phaseFuncFlags[[All,1]]]]; 
                     $Failed),
                Not @ MemberQ[bitEncodingFlags[[All,1]], OptionValue[BitEncoding]],
                    (Message[ApplyPhaseFunc::error, "Invalid option for BitEncoding. Must be one of " <> ToString@bitEncodingFlags[[All, 1]] <> "."]; 
                    $Failed),
                Not[ (overs === {}) || And[
                        MatchQ[overs, {({__Integer} -> _?Internal`RealValuedNumericQ) ..}],
                        Equal @@ Length /@ overs[[All,1]], 
                        Length[regs] === Length@overs[[1,1]] ] ],
                    (Message[ApplyPhaseFunc::error, "Invalid PhaseOverrides. Each overriden phase index must be specified as an n-tuple, where n is the number of sub-registers, pointing to a real number. For example, ApplyPhaseFunc[..., {{1},{2}}, ..., PhaseOverrides -> { {0,0} -> PI, ... }]."];
                     $Failed),
                StringEndsQ[func, "Distance"] && OddQ @ Length @ regs,
                    (Message[ApplyPhaseFunc::error, "'Distance' based phase functions require a strictly even number of subregisters, since every pair is assumed to represent the same coordinate."]; 
                    $Failed),
                Length[{params}] === 0,
                    ApplyNamedPhaseFuncInternal[qureg, Flatten[regs], Length/@regs, OptionValue[BitEncoding] /. bitEncodingFlags, func /. phaseFuncFlags, Flatten[overs[[All,1]]], N @ overs[[All,2]]],
                Length[{params}] > 0,
                    ApplyParamNamedPhaseFuncInternal[qureg, Flatten[regs], Length/@regs, OptionValue[BitEncoding] /. bitEncodingFlags, func /. phaseFuncFlags, N @ {params}, Flatten[overs[[All,1]]], N @ overs[[All,2]]]]]
        
        (* non-parameterised named func (multi-variable) *)
        ApplyPhaseFunc[qureg_Integer, regs:{{__Integer}..}, func_String, opts:OptionsPattern[]] := 
            ApplyPhaseFunc[qureg, regs, {func}, opts]
        
        (* invalid args and symbol syntax highlighting *)
        ApplyPhaseFunc[___] := invalidArgError[ApplyPhaseFunc]
        SyntaxInformation[ApplyPhaseFunc] = {"LocalVariables" -> {"Solve", {4, 4}}};
        
        
        
        (* 
         * Below is only front-end code for analytically simplifying expressions 
         * of Pauli tensors
         *)
         
         (* post-processing step to combine Pauli products that have identical symbols and indices... *)
        getPauliSig[ a: Subscript[(X|Y|Z|Id), _Integer] ] := {a}
        getPauliSig[ Verbatim[Times][t__] ] := Cases[{t}, Subscript[(X|Y|Z|Id), _]]
        getPauliSig[ _ ] := {}
        (* which works by splitting a sum into groups containing the same Pauli tensor, and simplifying each *)
        factorPaulis[s_Plus] := Total[Simplify /@ Plus @@@ GatherBy[List @@ s, getPauliSig]] /. Complex[0.`, 0.`] -> 0
        factorPaulis[e_] := e
        
        (* SimplifyPaulis prevents Mathemtica commutation (and inadvertently, variable substitution)
         * in order to perform all operator simplification correctly. It accepts expressions containing 
         * sums, products, powers and non-commuting multiples of Pauli operators (literals, and in variables), 
         * and the structures of remaining patterns/symbols/expressions can be anything at all.
         *)
        SetAttributes[SimplifyPaulis, HoldAll]
        
        SimplifyPaulis[ a:Subscript[(X|Y|Z|Id), _] ] := 
            a

        SimplifyPaulis[ (a:Subscript[(X|Y|Z), q_])^n_Integer ] /; (n >= 0) :=
        	If[EvenQ[n], Subscript[Id,q], a]
            
        SimplifyPaulis[ (a:Subscript[Id, q_])^n_Integer ] /; (n >= 0) :=
            a
        
        SimplifyPaulis[ Verbatim[Times][t__] ] := With[
        	(* hold each t (no substitutions), simplify each, release hold, simplify each (with subs) *)
        	{nc = SimplifyPaulis /@ (NonCommutativeMultiply @@ (SimplifyPaulis /@ Hold[t]))},
        	(* pass product (which now contains no powers of pauli expressions) to simplify *)
        	SimplifyPaulis[nc]]

        SimplifyPaulis[ Power[b_, n_Integer] ] /; (Not[FreeQ[b,Subscript[(X|Y|Z|Id), _]]] && n >= 0) :=
        	(* simplify the base, then pass a (non-expanded) product to simplify (to trigger above def) *)
        	With[{s=ConstantArray[SimplifyPaulis[b], n]}, 
        		SimplifyPaulis @@ (Times @@@ Hold[s])]
        		
        SimplifyPaulis[ Plus[t__] ] := With[
        	(* hold each t (no substitutions), simplify each, release hold, simplify each (with subs) *)
        	{s = Plus @@ (SimplifyPaulis /@ (List @@ (SimplifyPaulis /@ Hold[t])))},
        	(* combine identical Pauli tensors in the resulting sum *)
        	factorPaulis[s]
        ]

        SimplifyPaulis[ NonCommutativeMultiply[t__] ] := With[
        	(* hold each t (no substitutions), simplify each, release hold, simplify each (with subs) *)
        	{s = SimplifyPaulis /@ (NonCommutativeMultiply @@ (SimplifyPaulis /@ Hold[t]))},
        	(* expand all multiplication into non-commuting; this means ex can be a sum now *)
        	{ex = Distribute[s /. Times -> NonCommutativeMultiply]},
        	(* notation shortcuts *)
        	{xyz = X|Y|Z, xyzi = X|Y|Z|Id, ncm = NonCommutativeMultiply}, 
        	(* since ex can now be a sum, after below transformation, factorise *)
        	factorPaulis[
        		ex //. {
        		(* destroy exponents of single terms *)
        		(a:Subscript[xyzi, q_])^n_Integer  /; (n >= 0) :> If[EvenQ[n], Subscript[Id,q], a], 
                (a:Subscript[Id, _])^n_Integer :> a,
        		(* move scalars to their own element (to clean pauli pattern) *)
        		ncm[r1___, (f:Except[Subscript[xyzi, _]]) (a:Subscript[xyzi, _]) , r2___] :> ncm[f,r1,a,r2],
        		(* map same-qubit adjacent (closest) pauli matrices to their product *)
        		ncm[r1___, Subscript[(a:xyz), q_],r2:Shortest[___],Subscript[(b:xyz), q_], r3___] :>
        			If[a === b, ncm[r1,r2,r3,Subscript[Id,q]], With[{c = First @ Complement[List@@xyz, {a,b}]},
        				ncm[r1, I If[Sort[{a,b}]==={a,b},1,-1] Subscript[c, q], r2, r3]]],
                (* remove superfluous Id's when multiplying onto other paulis on any qubits *)
                ncm[r1___, Subscript[Id, _], r2___, b:Subscript[xyzi, _], r3___] :>
                    ncm[r1, r2, b, r3],
                ncm[r1___, a:Subscript[xyzi, _], r2___, Subscript[Id, _], r3___] :>
                    ncm[r1, a, r2, r3]
        	(* finally, restore products (overwriting user non-comms) and simplify scalars *)
        	} /. NonCommutativeMultiply -> Times]]
        	
        SimplifyPaulis[uneval_] := With[
            (* let everything else evaluate (to admit scalars, or let variables be substituted *)
            {eval=uneval},
            (* and if changed by its evaluation, attempt to simplify the new form *)
            If[Unevaluated[uneval] === eval, eval, SimplifyPaulis[eval]]]
            (* your eyes don't deceive you; === is smart enough to check inside Unevaluated! Nice one Stephen! *)
        
        
        
        (*
         * Below is only front-end code for generating 3D plots of density matrices
         *)
         
        Options[PlotDensityMatrix] = {
            BarSpacing -> .5,
            PlotComponent -> "Magnitude",
            ChartElementFunction -> "ProfileCube"
        };
         
        isNumericSquareMatrix[matrix_?MatrixQ] :=
            And[SquareMatrixQ @ matrix, AllTrue[Flatten @ matrix, NumericQ]]
        isNumericSquareMatrix[_] :=
            False
        isNumericVector[vector_?VectorQ] :=
            AllTrue[Flatten @ vector, NumericQ]
        isNumericVector[_] :=
            False
            
        extractMatrixData[comp_, matrix_] := With[
            {data = Switch[comp,
                "Real", Re @ matrix,
                "Imaginary", Im @ matrix,
                "Magnitude", Abs @ matrix,
                _, $Failed]},
            If[data === $Failed, Message[PlotDensityMatrix::error, "PlotComponent must be \"Real\", \"Imaginary\" or \"Magnitude\"."]];
            data]
        extractWeightedData[data_] :=
            WeightedData[Join @@ Array[List, Dimensions@#], Join @@ #]& @ Abs @ data (* forced positive *)
        extractChartElemFunc[func_] := With[
            {elem=(func /. Automatic -> 
                OptionValue[Options[PlotDensityMatrix], ChartElementFunction])},
            If[Not @ StringQ @ elem, 
                Message[PlotDensityMatrix::error, "ChartElementFunction must be a string, or Automatic. See available options with ChartElementData[Histogram3D]."]; 
                $Failed,
                elem]]
        extractBarSpacing[val_] := With[
            {space = (val /. Automatic -> OptionValue[Options[PlotDensityMatrix], BarSpacing])},
            If[Not[NumericQ@space] || Not[0 <= space < 1],
                Message[PlotDensityMatrix::error, "BarSpacing (Automatic is .5) must be a number between 0 (inclusive) and 1 (exclusive)."];
                $Failed,
                space]]
            
        plotDensOptsPatt = OptionsPattern[{PlotDensityMatrix,Histogram3D}];
         
        (* single matrix plot *)
        PlotDensityMatrix[id_Integer, opts:plotDensOptsPatt] :=
            PlotDensityMatrix[GetQuregMatrix[id], opts]
        PlotDensityMatrix[matrix_?isNumericSquareMatrix, opts:plotDensOptsPatt] :=
            Block[{data, chartelem, space, offset},
                (* unpack data and args (may throw Message) *)
                data = extractMatrixData[OptionValue[PlotComponent], matrix];
                chartelem = extractChartElemFunc[OptionValue[ChartElementFunction]];
                space = extractBarSpacing[OptionValue[BarSpacing]];
                offset = space {1,-1}/2;
                (* return early if error *)
                If[MemberQ[{data,chartelem,space}, $Failed], Return[$Failed]];
                (* plot *)
                Histogram3D[
                    extractWeightedData[data], 
                    Times @@ Dimensions @ data,
                    (* offset and possibly under-axis (negative values) bar graphics *)
                    ChartElementFunction->
            			Function[{region, inds, meta}, ChartElementData[chartelem][
                            (* we subtract "twice" negative values, to point e.g. cones downard *)
                            region + {offset, offset, {0, 2 Min[0, Extract[data, First@inds]]}}, inds, meta]],
                    (* with user Histogram3D options *)
                    FilterRules[{opts}, Options[Histogram3D]],
                    (* and our overridable defaults *)
                    ColorFunction -> (ColorData["DeepSeaColors"][1 - #] &),
                    PlotRange -> {
                        .5 + {0, First @ Dimensions @ data},
                        .5 + {0, Last @ Dimensions @ data},
                        Automatic}
                ]
            ]
        (* two matrix plot *)
        PlotDensityMatrix[id1_Integer, id2_Integer, opts:plotDensOptsPatt] :=
            PlotDensityMatrix[GetQuregMatrix[id1], GetQuregMatrix[id2], opts]
        PlotDensityMatrix[id1_Integer, matr2_?isNumericSquareMatrix, opts:plotDensOptsPatt] :=
            PlotDensityMatrix[GetQuregMatrix[id1], matr2, opts]
        PlotDensityMatrix[matr1_?isNumericSquareMatrix, id2_Integer, opts:plotDensOptsPatt] :=
            PlotDensityMatrix[matr1, GetQuregMatrix[id2], opts]
        PlotDensityMatrix[matr1_?isNumericSquareMatrix, vec2_?isNumericVector, opts:plotDensOptsPatt] :=
            With[{matr2 = KroneckerProduct[ConjugateTranspose@{vec2}, vec2]},
                PlotDensityMatrix[matr1, matr2, opts]]
        PlotDensityMatrix[matr1_?isNumericSquareMatrix, matr2_?isNumericSquareMatrix, opts:plotDensOptsPatt] :=
            Block[{data1, data2, chartelem, space, offset},
                (* unpack data and args (may throw Message) *)
                data1 = extractMatrixData[OptionValue[PlotComponent], matr1];
                data2 = extractMatrixData[OptionValue[PlotComponent], matr2];
                chartelem = extractChartElemFunc[OptionValue[ChartElementFunction]];
                space = extractBarSpacing[OptionValue[BarSpacing]];
                offset = space {1,-1}/2;
                (* return early if error *)
                If[MemberQ[{data1,data2,chartelem,space}, $Failed], Return[$Failed]];
                (* plot *)
                Histogram3D[
                    {extractWeightedData[data1], extractWeightedData[data2]}, 
                    Max[Times @@ Dimensions @ data1, Times @@ Dimensions @ data2],
                    (* offset and possibly under-axis (negative values) bar graphics *)
                    ChartElementFunction -> {
                        Function[{region, inds, meta}, ChartElementData[chartelem][
                            (* we subtract "twice" negative values, to point e.g. cones downard *)
                            region + {offset, offset, {0, 2 Min[0, Extract[data1, First@inds]]}}, inds, meta]],
                        Function[{region, inds, meta}, ChartElementData[chartelem][
                            (* we make the second matrix bars slightly inside the first's *)
                            region + {1.001 offset, 1.001 offset, {0, 2 Min[0, Extract[data2, First@inds]]}}, inds, meta]]
                        },
                    (* with user Histogram3D options *)
                    FilterRules[{opts}, Options[Histogram3D]],
                    (* and our overridable defaults *)
                    ChartStyle -> {Opacity[1], Opacity[.3]},
                    ColorFunction -> (ColorData["DeepSeaColors"][1 - #] &),
                    PlotRange -> {
                        .5 + {0, Max[First @ Dimensions @ data1, First @ Dimensions @ data2]},
                        .5 + {0, Max[Last @ Dimensions @ data1, Last @ Dimensions @ data2]},
                        Automatic},
                    (* useless placebo *)
                    Method -> {"RelieveDPZFighting" -> True}
                ]
            ]
        PlotDensityMatrix[___] := (
            Message[PlotDensityMatrix::error, "Invalid arguments. See ?PlotDensityMatrix. Note the first argument must be a numeric square matrix."]; 
            $Failed)
        
        
        
        (*
         * Below is only front-end code for generating circuit diagrams from
         * from circuit the same format circuit specification
         *)
         
        (* convert symbolic gate form to {symbol, ctrls, targets} *)
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][ R[arg_, Verbatim[Times][paulis:Subscript[pauliCodePatt, _Integer]..]] ]] := {Join[{R}, {paulis}[[All,1]]], {ctrls}, {paulis}[[All,2]]}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][ R[arg_, Subscript[pauli:pauliCodePatt, targ_Integer]] ]] := {{R,pauli}, {ctrls}, {targ}}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]] := {gate, {},{targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]] := {gate, {}, {targs}}
        getSymbCtrlsTargs[R[arg_, Verbatim[Times][paulis:Subscript[pauliCodePatt, _Integer]..]]] := {Join[{R}, {paulis}[[All,1]]], {}, {paulis}[[All,2]]}
        getSymbCtrlsTargs[R[arg_, Subscript[pauli:(X|Y|Z), targ_Integer]]] := {{R,pauli}, {}, {targ}}
            (* little hack to enable G[x] and Fac[y] in GetCircuitColumns *)
            getSymbCtrlsTargs[G[x_]] := {G, {}, {}}
            getSymbCtrlsTargs[Fac[x_]] := {Fac, {}, {}}

        (* deciding how to handle gate placement *)
        getQubitInterval[{ctrls___}, {targs___}] :=
        	Interval @ {Min[ctrls,targs],Max[ctrls,targs]}
        getNumQubitsInCircuit[circ_List] :=
        	Max[1 + Cases[{circ}, Subscript[gate_, inds__]-> Max[inds], Infinity],    
        		1 + Cases[{circ}, Subscript[gate_, inds__][___] -> Max[inds], Infinity]] /. -Infinity -> 1  (* assume G and Fac circuits are 1 qubit *)
        isContiguousBlockGate[(SWAP|M|Rz|Ph|X|R|{R, pauliCodePatt..})] := False
        isContiguousBlockGate[_] := True
        needsSpecialSwap[label_, _List] /; Not[isContiguousBlockGate[label]] := False
        needsSpecialSwap[label_Symbol, targs_List] :=
        	And[Length[targs] === 2, Abs[targs[[1]] - targs[[2]]] > 1]
        getFixedThenBotTopSwappedQubits[{targ1_,targ2_}] :=
        	{Min[targ1,targ2],Min[targ1,targ2]+1,Max[targ1,targ2]}
            
        (* gate and qubit graphics primitives *)
        drawCross[targ_, col_] := {
        	Line[{{col+.5,targ+.5}-{.1,.1},{col+.5,targ+.5}+{.1,.1}}],
        	Line[{{col+.5,targ+.5}-{-.1,.1},{col+.5,targ+.5}+{-.1,.1}}]}
        drawControls[{ctrls__}, {targs___}, col_] := {
        	FaceForm[Black],
        	Table[Disk[{col+.5,ctrl+.5},.1],{ctrl,{ctrls}}],
        	With[{top=Max@{ctrls,targs},bot=Min@{ctrls,targs}},
        		Line[{{col+.5,bot+.5},{col+.5,top+.5}}]]}
        drawSingleBox[targ_, col_] :=
        	Rectangle[{col+.1,targ+.1}, {col+1-.1,targ+1-.1}]
        drawDoubleBox[targ_, col_] :=
        	Rectangle[{col+.1,targ+.1}, {col+1-.1,targ+2-.1}]
        drawMultiBox[minTarg_, numTargs_, col_] :=
            Rectangle[{col+.1,minTarg+.1}, {col+1-.1,minTarg+numTargs-.1}]
        drawQubitLines[qubits_List, col_, width_:1] :=
        	Table[Line[{{col,qb+.5},{col+width,qb+.5}}], {qb,qubits}]
        drawSpecialSwapLine[targ1_, targ2_, col_] := {
        	Line[{{col,targ1+.5},{col+.1,targ1+.5}}],
        	Line[{{col+.1,targ1+.5},{col+.5-.1,targ2+.5}}],
        	Line[{{col+.5-.1,targ2+.5},{col+.5,targ2+.5}}]}
        drawSpecialSwap[targ1_,targ2_,col_] := {
        	drawSpecialSwapLine[targ1,targ2,col],
        	drawSpecialSwapLine[targ2,targ1,col]}
            
        (* special gate graphics *)
        drawGate[SWAP, {}, {targs___}, col_] := {
            (drawCross[#,col]&) /@ {targs},
            Line[{{col+.5,.5+Min@targs},{col+.5,.5+Max@targs}}]}
        drawGate[Z, {ctrls__}, {targ_}, col_] := {
            drawControls[{ctrls,targ},{targ},col],
            Line[{{col+.5,.5+Min@ctrls},{col+.5,.5+Max@ctrls}}]}
        drawGate[Ph, {ctrls___}, {targs__}, col_] := {
            drawControls[{ctrls,targs},{},col],
            Text["\[Theta]", {col+.75,Min[{ctrls,targs}]+.75}]}
        drawGate[label:(Kraus|KrausNonTP|Damp|Deph|Depol), {}, targs_List, col_] := {
            EdgeForm[Dashed],
            drawGate[label /. {
                    Kraus -> \[Kappa], KrausNonTP -> \[Kappa]NTP, Damp -> \[Gamma], 
                    Deph -> \[Phi], Depol -> \[CapitalDelta]},
                {}, targs, col]}

        (* single qubit gate graphics *)
        drawGate[Id, {}, {targs___}, col_] :=
            {}
        drawGate[visibleId, {}, {targ_}, col_] := {
            drawSingleBox[targ, col],
            Text["\[DoubleStruckOne]", {col+.5,targ+.5}]
        }
        drawGate[visibleId, {}, {targs___}, col_] :=
            drawGate[visibleId, {}, #, col]& /@ {targs}
        drawGate[M, {}, {targs___}, col_] :=
        	Table[{
        		drawSingleBox[targ,col],
        		Circle[{col+.5,targ+.5-.4}, .4, {.7,\[Pi]-.7}],
        		Line[{{col+.5,targ+.5-.25}, {col+.5+.2,targ+.5+.3}}]
        		}, {targ, {targs}}]

        drawGate[Depol, {}, {targ_}, col_] := {
            EdgeForm[Dashed], drawGate[\[CapitalDelta], {}, {targ}, col]}
        drawGate[X, {}, {targ_}, col_] := {
            Circle[{col+.5,targ+.5},.25],
            Line[{{col+.5,targ+.5-.25},{col+.5,targ+.5+.25}}]}
        drawGate[label_Symbol, {}, {targ_}, col_] := {
        	drawSingleBox[targ, col],
        	Text[SymbolName@label, {col+.5,targ+.5}]}
            
        (* multi-qubit gate graphics *)
        drawGate[Rz, {}, targs_List, col_] := {
            Line[{{col+.5,Min[targs]+.5},{col+.5,Max[targs]+.5}}],
            Sequence @@ (drawGate[Rz, {}, {#1}, col]& /@ targs)}
        drawGate[{R, rots:pauliCodePatt..}, {}, targs_List, col_] := {
            Line[{{col+.5,Min[targs]+.5},{col+.5,Max[targs]+.5}}],
            Sequence @@ MapThread[drawGate[#1/.{X->Rx,Y->Ry,Z->Rz,Id->visibleId}, {}, {#2}, col]&, {{rots}, targs}]}
        drawGate[G, {}, targs_List, col_] /; (isContiguousBlockGate[label] && Union@Differences@Sort@targs=={1}) := {
            drawMultiBox[Min[targs], Length[targs], col],
            Text["e"^"i\[Theta]", {col+.5,Mean[targs]+.5}]}
        drawGate[Fac, {}, targs_List, col_] /; (isContiguousBlockGate[label] && Union@Differences@Sort@targs=={1}) := {
            drawMultiBox[Min[targs], Length[targs], col],
            Text[Rotate["factor",Pi/2], {col+.5,Mean[targs]+.5}]}
        drawGate[label_Symbol, {}, targs_List, col_] /; (isContiguousBlockGate[label] && Union@Differences@Sort@targs=={1}) := {
            drawMultiBox[Min[targs], Length[targs], col],
            Text[SymbolName@label, {col+.5,Mean[targs]+.5}]}
        drawGate[label_Symbol, {}, targs_List, col_] := {
            Line[{{col+.5,Min[targs]+.5},{col+.5,Max[targs]+.5}}],
            Sequence @@ (drawGate[label, {}, {#1}, col]& /@ targs)}
                
        (* two-qubit gate graphics *)
        drawGate[X, {}, targs:{targ1_,targ2_}, col_] := {
            Line[{{col+.5,targ1+.5},{col+.5,targ2+.5}}],
            Sequence @@ (drawGate[X, {}, {#1}, col]& /@ targs)}
        drawGate[label_Symbol, {}, {targ1_,targ2_}/;Abs[targ2-targ1]===1, col_] := {
        	drawDoubleBox[Min[targ1,targ2], col],
        	Text[SymbolName@label, {col+.5,Min[targ1,targ2]+.5+.5}]}
        drawGate[label_Symbol, {}, {targ1_,targ2_}, col_] := 
        	With[{qb = getFixedThenBotTopSwappedQubits[{targ1,targ2}]}, {
        		drawSpecialSwap[qb[[3]], qb[[2]], col],
        		drawGate[label, {}, {qb[[1]],qb[[2]]}, col+.5],
        		drawSpecialSwap[qb[[3]],qb[[2]],col+.5+1]}]
        		
        (* controlled gate graphics *)
        drawGate[SWAP, {ctrls__}, {targs__}, col_] := {
        	drawControls[{ctrls},{targs},col],
        	drawGate[SWAP, {}, {targs}, col]}
        drawGate[ops:{R, (X|Y|Z)..}, {ctrls__}, {targs__}, col_] := {
            drawControls[{ctrls},{targs},col],
            drawGate[ops, {}, {targs}, col]}
        drawGate[label_Symbol, {ctrls__}, {targ_}, col_] := {
        	drawControls[{ctrls},{targ},col],
        	drawGate[label, {}, {targ}, col]}
        drawGate[label_Symbol, {ctrls__}, targs_List, col_] := {
        	If[needsSpecialSwap[label,targs],
        		With[{qb = getFixedThenBotTopSwappedQubits[targs]},
        			drawControls[{ctrls} /. qb[[2]]->qb[[3]], DeleteCases[targs,qb[[3]]], col+.5]],
        		drawControls[{ctrls},targs, col]],
        	drawGate[label, {}, targs, col]}
        
        (* generating background qubit lines in a whole circuit column *)
        drawQubitColumn[isSpecialSwapCol_, specialSwapQubits_, numQubits_, curCol_, width_:1] :=
        	If[isSpecialSwapCol,
        		(* for a special column, draw all middle lines then non-special left/right nubs *)
        		With[{nonspecial=Complement[Range[0,numQubits-1], specialSwapQubits]}, {
        			drawQubitLines[Range[0,numQubits-1],curCol+.5,width],
        			drawQubitLines[nonspecial,curCol,width],
        			drawQubitLines[nonspecial,curCol+1,width]}],
        		(* for a non special column, draw all qubit lines *)
        		drawQubitLines[Range[0,numQubits-1],curCol,width]
        	]
            
        (* subcircuit seperator and scheduling label graphics *)
        drawSubcircSpacerLine[col_, numQubits_, style_] := 
            If[style === None, Nothing, {
                style,
                Line[{{col,0},{col,numQubits}}]}]
            
        (* labels below seperator *)
        defaultSubcircLabelDrawer[label_, col_] :=
            Rotate[Text[label, {col,-.5}], -15 Degree]
        defaultTimeLabelDrawer[time_, col_] := 
            defaultSubcircLabelDrawer[DisplayForm @ RowBox[{"t=",time}], col]
                
        (* optionally compactifies a circuit via GetColumnCircuits[] *)
        compactCirc[flag_][circ_] :=
            If[flag, Flatten @ GetCircuitColumns[circ], circ]
            
        (* creates a graphics description of the given list of subcircuits *)
        generateCircuitGraphics[subcircs_List, numQubits_Integer, opts___] := With[{
            
            (* unpack optional arguments *)
            subSpacing=OptionValue[{opts,Options[DrawCircuit]}, SubcircuitSpacing],
            dividerStyle=OptionValue[{opts,Options[DrawCircuit]}, DividerStyle],
            labelDrawFunc=OptionValue[{opts,Options[DrawCircuit]}, LabelDrawer],
            labels=OptionValue[{opts, Options[DrawCircuit]}, SubcircuitLabels],
            compactFlag=OptionValue[{opts,Options[DrawCircuit]}, Compactify]
            },
            
            (* prepare local variables *)
            Module[{
            	qubitgraphics,gategraphics,
            	curCol,curSymb,curCtrls,curTargs,curInterval,curIsSpecialSwap,
            	prevIntervals,prevIsSpecialSwap,prevSpecialQubits,
                isFirstGate,subcircInd=0,subcircIndMax=Length[subcircs],
                gates},

                (* outputs *)
            	qubitgraphics = {};
            	gategraphics = {};
            	
                (* status of whether a gate can fit in the previous column *)
            	curCol = 0;
            	prevIntervals = {};
            	prevIsSpecialSwap = False;
            	prevSpecialQubits = {};
                
                (* draw divider left-of-circuit for the first label (else don't draw it) *)
                If[And[labelDrawFunc =!= None, labels =!= {}, First[labels] =!= None], 
                    AppendTo[gategraphics, drawSubcircSpacerLine[subSpacing/2, numQubits, dividerStyle]];
                    AppendTo[gategraphics, labelDrawFunc[First @ labels, subSpacing/2]];
                    AppendTo[qubitgraphics, drawQubitColumn[False, {} , numQubits, curCol + subSpacing/2, subSpacing/2]];
                    curCol = subSpacing;
                ];
                
                (* for each subcircuit... *)
                Table[
                    gates = compactCirc[compactFlag][subcirc /. {
                        (* hackily replace qubit-free gates full-state contiguous gate *)
                        G[_] -> Subscript[G, Range[0,numQubits-1] ], 
                        Fac[_] -> Subscript[Fac, Range[0,numQubits-1] ]}];
                    subcircInd++;
                    isFirstGate = True;
                
                    (* for each gate... *)
                	Table[
                    
                		(* extract data from gate *)
                		{curSymb,curCtrls,curTargs} = getSymbCtrlsTargs[gate];
                		curInterval = getQubitInterval[curCtrls,curTargs];
                		curIsSpecialSwap = needsSpecialSwap[curSymb,curTargs];
                		
                		(* decide whether gate will fit in previous column *)
                		If[Or[
                			And[curIsSpecialSwap, Not[prevIsSpecialSwap]],
                			AnyTrue[prevIntervals, (IntervalIntersection[curInterval,#] =!= Interval[]&)]],
                			(* will not fit: *)
                			(

                				(* draw qubit lines for the previous column (if it exists) *)
                                If[Not[isFirstGate],
                                    AppendTo[qubitgraphics,
                                        drawQubitColumn[prevIsSpecialSwap, prevSpecialQubits, numQubits, curCol]]];
                				
                				(* then make a new column *)
                				curCol = curCol + If[isFirstGate, 0, If[prevIsSpecialSwap,2,1]]; 
                				prevIntervals = {curInterval};
                				prevIsSpecialSwap = curIsSpecialSwap;
                				prevSpecialQubits = {};
                                    
                			),
                			(* will fit *)
                			AppendTo[prevIntervals, curInterval]
                		];
                        
                        (* record that this is no longer the first gate in the subcircuit *)
                        isFirstGate = False;
                		
                		(* record whether this gate requires special swaps *)
                		If[curIsSpecialSwap, 
                			With[{qbs=getFixedThenBotTopSwappedQubits[curTargs]},
                				AppendTo[prevSpecialQubits, qbs[[2]]];
                				AppendTo[prevSpecialQubits, qbs[[3]]]]];
                	
                		(* draw gate *)
                		AppendTo[gategraphics,
                			drawGate[
                				curSymb,curCtrls,curTargs,
                				curCol + If[prevIsSpecialSwap ~And~ Not[curIsSpecialSwap], .5, 0]]];        
                        ,
                		{gate,gates}
                	];
            	
                	(* perform the final round of qubit line drawing (for previous column) *)
                	AppendTo[qubitgraphics, 
                		drawQubitColumn[prevIsSpecialSwap, prevSpecialQubits, numQubits, curCol]];
                        
                    (* make a new column (just accounting for previous subcircuit) *)
                    curCol = curCol + If[prevIsSpecialSwap,2,1]; 
                    prevIntervals = {};
                    prevIsSpecialSwap = False;
                    prevSpecialQubits = {};
                        
                    (* if this was not the final subcircuit... *)
                    If[subcircInd < subcircIndMax, 
                    
                        (* add subcircit seperator line  *)
                        AppendTo[gategraphics, 
                            drawSubcircSpacerLine[curCol + subSpacing/2, numQubits, dividerStyle]];
                                
                        (* add label below seperator line (unless None) *)
                        If[And[labelDrawFunc =!= None, subcircInd+1 <= Length[labels], labels[[subcircInd+1]] =!= None],
                            AppendTo[gategraphics, 
                                labelDrawFunc[labels[[subcircInd+1]], curCol + subSpacing/2]]];
                    
                        (* add offset for subcircuit spacing (avoid 0 length for visual artefact) *)
                        If[subSpacing > 0,
                            AppendTo[qubitgraphics, 
                                drawQubitColumn[prevIsSpecialSwap, prevSpecialQubits, numQubits, curCol, subSpacing]];
                            curCol = curCol + subSpacing];
                    ]
                    ,
                    {subcirc, subcircs}
                ];
                
                (* if there's a remaining label, draw it and a final (otherwise superfluous) divider *)
                If[And[labelDrawFunc =!= None, subcircInd+1 <= Length[labels], labels[[subcircInd+1]] =!= None],
                    AppendTo[gategraphics, 
                        labelDrawFunc[labels[[subcircInd+1]], curCol + subSpacing/2]];
                    AppendTo[gategraphics, 
                        drawSubcircSpacerLine[curCol + subSpacing/2, numQubits, dividerStyle]];
                    AppendTo[qubitgraphics,
                        drawQubitColumn[False, {} , numQubits, curCol, subSpacing/2]];
                ];
            	
                (* return *)
            	{curCol, qubitgraphics, gategraphics}
            ]    
        ]
        
        (* renders a circuit graphics description *)
        displayCircuitGraphics[{numCols_, qubitgraphics_, gategraphics_}, opts___] :=
            Graphics[
                {
                    FaceForm[White], EdgeForm[Black],
                    qubitgraphics, gategraphics
                },
                FilterRules[{opts}, Options[Graphics]],
                ImageSize -> 30 (numCols+1),
                PlotRangePadding -> None
                (* , Method -> {"ShrinkWrap" -> True} *)
            ]
        
        (* declaring optional args to DrawCircuit *)
        Options[DrawCircuit] = {
            Compactify -> True,
            DividerStyle -> Directive[Dashed, Gray], (* None for no dividers *)
            SubcircuitSpacing -> .25,
            SubcircuitLabels -> {},
            LabelDrawer -> defaultSubcircLabelDrawer
        };
        
        (* public functions to fully render a circuit *)
        DrawCircuit[noisySched:{{_, _List, _List}..}, Repeated[numQubits_, {0,1}], opts:OptionsPattern[{DrawCircuit,Graphics}]] :=
            (* compactify each subcirc but not their union *)
            DrawCircuit[{First[#], Join @@ compactCirc[OptionValue[Compactify]] /@ Rest[#]}& /@ noisySched, numQubits, Compactify -> False, opts]
        DrawCircuit[schedule:{{_, _List}..}, numQubits_Integer, opts:OptionsPattern[]] :=
            displayCircuitGraphics[
                generateCircuitGraphics[schedule[[All,2]], numQubits, opts, SubcircuitLabels -> schedule[[All,1]], LabelDrawer -> defaultTimeLabelDrawer], opts]
        DrawCircuit[schedule:{{_, _List}..}, opts:OptionsPattern[]] :=
            DrawCircuit[schedule, getNumQubitsInCircuit[Flatten @ schedule[[All,2]]], opts]
        DrawCircuit[cols:{{___}..}, numQubits_Integer, opts:OptionsPattern[]] := 
            displayCircuitGraphics[
                generateCircuitGraphics[cols, numQubits, {}, opts], opts]
        DrawCircuit[cols:{{___}..}, opts:OptionsPattern[]] :=
            DrawCircuit[cols, getNumQubitsInCircuit[Flatten @ cols], opts]
        DrawCircuit[circ_List, args___] :=
            DrawCircuit[{circ}, args]
        DrawCircuit[___] := invalidArgError[DrawCircuit]
        
        
        
        (*
         * Below is only front-end code for drawing topology diagrams of circuits
         *)
         
        (* every gate is distinguished in 'Parameters' mode, even if differing only by parameter *)
        getTopolGateLabel["Parameters"][g_] := g
        
        (* in 'Qubits' mode, parameters are removed/ignored *)
        getTopolGateLabel["Qubits"][Subscript[C, q__][s_]] := Subscript[C, q][ getTopolGateLabel["Qubits"][s] ]
        getTopolGateLabel["Qubits"][R[_, p_]] := R[p]
        getTopolGateLabel["Qubits"][Subscript[s_, q__][__]] := Subscript[s, q]
        getTopolGateLabel["Qubits"][s_] := s
        
        (* in 'NumQubits' mode, specific qubit indices are ignored, replaced by a-z *)
        getLetterSeq[len_Integer, start_] := Sequence @@ Table[ FromLetterNumber[t], {t,start,start+len-1} ]
        getLetterSeq[list_List, start_] := getLetterSeq[Length[list], start]
        alphabetisizeGate[Subscript[C, {q__}|q__][s_]] := Subscript[C, getLetterSeq[{q},1]][ alphabetisizeGate[s, 1+Length@{q}]]
        alphabetisizeGate[Subscript[s_, {q__}|q__], i_:1] := Subscript[s, getLetterSeq[{q},i]]
        alphabetisizeGate[Subscript[s_, {q__}|q__][_], i_:1] := Subscript[s, getLetterSeq[{q},i]]
        alphabetisizeGate[R[_, p_Times], i_:1] := R[Times @@ Subscript @@@ Transpose[{(List @@ p)[[All,1]], List @ getLetterSeq[Length[p],i]}]]
        alphabetisizeGate[R[_, Subscript[p_, q_]], i_:1] := R[Subscript[p, getLetterSeq[1,i]]]
        getTopolGateLabel["NumberOfQubits"][g_] := alphabetisizeGate[g]
        
        (* in 'Gates' mode, all qubits are ignored and removed *)
        getTopolGateLabel["Gates"][Subscript[C, q__][s_]] := C[ getTopolGateLabel["Gates"][s] ]
        getTopolGateLabel["Gates"][R[_, p_Times]] := R[Row @ (List @@ p)[[All,1]]]
        getTopolGateLabel["Gates"][R[_, Subscript[p_, q_]]] := R[p]
        getTopolGateLabel["Gates"][Subscript[s_, __]|Subscript[s_, __][__]] := s
        
        (* in 'None' mode, all gate properties are discarded *)
        getTopolGateLabel["None"][s_] := None
        
        (* 'Connectivity' mode is handled entirely in DrawCircuitTopology[] *)
        
        getCircTopolGraphData[circ_, showReps_, showLocal_, groupMode_] := Module[
            {edges = {}, edgeLabels = <||>},
            
            (* for every gate ... *)
            Table[ With[
                (* extract targeted qubits, choose group/label *)
                {qubits = getSymbCtrlsTargs[gate][[{2,3}]] // Flatten},
                {label = If[groupMode === "Connectivity", 
                    Row[Sort[qubits],Spacer[.1]], 
                    getTopolGateLabel[groupMode][gate]]},
                    
                (* optionally admit single-qubit gates *)
                {vertices = If[showLocal && Length[qubits] === 1, 
                    Join[qubits, qubits],
                    qubits]},
                
                (* for every pair of targeted qubits... *)
                Table[ With[ {key = UndirectedEdge @@ Sort[pair] },
                    If[ KeyExistsQ[edgeLabels, key],
                    
                        (* if the edge exists, but the label is new (or we allow repetition), record it (else do nothing) *)
                        If[ showReps || Not @ MemberQ[edgeLabels[key], label],
                            AppendTo[edges, key];
                            AppendTo[edgeLabels[key], label]],
                            
                        (* else if the edge is new, record it unconditionally *)
                        AppendTo[edges, key];
                        edgeLabels[key] = {label}
                    ]],
                    {pair, Subsets[vertices, {2}]}
                ]],
                {gate, circ}];
                
            (* return *)
            {edges, edgeLabels}]
        
        Options[DrawCircuitTopology] = {
            ShowRepetitions -> False,
            ShowLocalGates -> True,
            DistinguishBy -> "Gates",
            DistinguishedStyles -> Automatic
        };
        
        DrawCircuitTopology[circ_List, opts:OptionsPattern[{DrawCircuitTopology, Graph, Show, LineLegend}]] := Module[
            {edges, edgeLabels, edgeIndices, graph},
            
            (* validate opt args *)
            
            (* extract topology data from circuit *)
            {edges, edgeLabels} = getCircTopolGraphData[circ, 
                OptionValue[{opts,Options[DrawCircuitTopology]}, ShowRepetitions], 
                OptionValue[{opts,Options[DrawCircuitTopology]}, ShowLocalGates], 
                OptionValue[{opts,Options[DrawCircuitTopology]}, DistinguishBy]];
            
            (* maintain indices for the list of labels recorded for each edge (init to 1) *)
            edgeIndices = <| Rule @@@ Transpose[{Keys[edgeLabels], ConstantArray[1, Length[edgeLabels]]}] |>;
        
            (* prepare an undirected graph (graphical form) with place-holder styles *)
            graph = Show[ 
                Graph[
                    Table[Style[edge, STYLES[edge]], {edge,edges}],
                    (* user-overridable default Graph properties *)
                    Sequence @@ FilterRules[{opts}, Options[Graph]],     (* Sequence[] shouldn't be necessary; another MMA Graph bug, sigh! *)
                    VertexStyle -> White,
                    VertexSize -> .1,
                    VertexLabels -> Automatic],
                (* user-overridable Show options *)
                FilterRules[{opts}, Options[Show]]];
                
            (* if there are no distinguishing gate groups, remove place-holder styles and return graph *)
            If[OptionValue[DistinguishBy] === "None",
                Return[ graph /. STYLES[_] -> Automatic ]];
                
            (* otherwise...  *)
            With[
                (* prepare legend with unique labels *)
                {legLabels = DeleteDuplicates @ Flatten @ Values @ edgeLabels},
                {legStyles = If[
                        OptionValue[DistinguishedStyles] === Automatic,
                        ColorData["Rainbow"] /@ Range[0,1,1/(Length[legLabels]-1)],
                        PadRight[OptionValue[DistinguishedStyles], Length[legLabels], OptionValue[DistinguishedStyles]]
                ]},
                {edgeStyles = Rule @@@ Transpose[{legLabels, legStyles}]},
                
                (* style each graph edge one by one *)
                {graphic = graph /. STYLES[edge_] :> (edgeLabels[edge][[ edgeIndices[edge]++ ]] /. edgeStyles)},
                
                (* return the styled graph with a line legend *)
                Legended[graphic,
                    LineLegend[legStyles, StandardForm /@ legLabels, 
                        Sequence @@ FilterRules[{opts}, Options[LineLegend]]]]
            ]
        ]
                
        
        
        (*
         * Below are front-end functions for 
         * generating analytic expressions from
         * circuit specifications
         *)
         
        (* generate a swap matrix between any pair of qubits *)
        getAnalSwapMatrix[qb1_, qb1_, numQb_] :=
            IdentityMatrix @ Power[2,numQb]
        getAnalSwapMatrix[qb1_, qb2_, numQb_] /; (qb1 > qb2) :=
    	    getAnalSwapMatrix[qb2, qb1, numQb]
        getAnalSwapMatrix[qb1_, qb2_, numQb_] := Module[
           	{swap, iden, p0,l1,l0,p1, block=Power[2, qb2-qb1]},
           	
           	(* as per Lemma 3.1 of arxiv.org/pdf/1711.09765.pdf *)
           	iden = IdentityMatrix[block/2];
           	p0 = KroneckerProduct[iden, {{1,0},{0,0}}];
           	l0 = KroneckerProduct[iden, {{0,1},{0,0}}];
           	l1 = KroneckerProduct[iden, {{0,0},{1,0}}];
           	p1 = KroneckerProduct[iden, {{0,0},{0,1}}];
           	swap = ArrayFlatten[{{p0,l1},{l0,p1}}];
           	
           	(* pad swap matrix to full size *)
           	If[qb1 > 0, 
           		swap = KroneckerProduct[swap, IdentityMatrix@Power[2,qb1]]];
           	If[qb2 < numQb-1, 
           		swap = KroneckerProduct[IdentityMatrix@Power[2,numQb-qb2-1], swap]];
           	swap
        ]
        
        (* build a numQb-matrix from op matrix *)
        getAnalFullMatrix[ctrls_, targs_, op_, numQb_] := Module[
        	{swaps=IdentityMatrix@Power[2,numQb], unswaps, swap, newctrls, newtargs, i,j, matr},
            
            (* make copies of user lists so we don't modify *)
        	unswaps = swaps;
        	newctrls = ctrls;
        	newtargs = targs;
        	
        	(* swap targs to {0,1,...} *)
        	For[i=0, i<Length[newtargs], i++,
        		If[i != newtargs[[i+1]],
        			swap = getAnalSwapMatrix[i, newtargs[[i+1]], numQb];
        			swaps = swap . swaps;
        			unswaps = unswaps . swap;
        			newctrls = newctrls /. i->newtargs[[i+1]];
        			newtargs = newtargs /. i->newtargs[[i+1]];
        		]
        	];
        	
        	(* swap ctrls to {Length[targs], Length[targs]+1, ...} *)
        	For[i=0, i<Length[newctrls], i++,
        		j = Length[newtargs] + i;
        		If[j != newctrls[[i+1]],
        			swap = getAnalSwapMatrix[j, newctrls[[i+1]], numQb];
        			swaps = swap . swaps;
        			unswaps = unswaps . swap;
        			newctrls = newctrls /. j->newctrls[[i+1]]
        		]
        	];
        	
        	(* create controlled(op) *)
        	matr = IdentityMatrix @ Power[2, Length[ctrls]+Length[targs]];
        	matr[[-Length[op];;,-Length[op];;]] = op;
        	
        	(* pad to full size *)
        	If[numQb > Length[targs]+Length[ctrls],
        		matr = KroneckerProduct[
        			IdentityMatrix@Power[2, numQb-Length[ctrls]-Length[targs]], matr]
        	];
        	
        	(* apply the swaps on controlled(op) *)
        	matr = unswaps . matr . swaps;
        	matr
        ]
        
        (* map gate symbols to matrices *)
        getAnalGateMatrix[Subscript[H, _]] = {{1,1},{1,-1}}/Sqrt[2];
        getAnalGateMatrix[Subscript[X, _]] = PauliMatrix[1];
        getAnalGateMatrix[Subscript[X, t__]] := KroneckerProduct @@ ConstantArray[PauliMatrix[1],Length[{t}]];
        getAnalGateMatrix[Subscript[Y, _]] = PauliMatrix[2];
        getAnalGateMatrix[Subscript[Z, _]] = PauliMatrix[3];
        getAnalGateMatrix[Subscript[S, _]] = {{1,0},{0,I}};
        getAnalGateMatrix[Subscript[T, _]] = {{1,0},{0,Exp[I Pi/4]}};
        getAnalGateMatrix[Subscript[SWAP, _,_]] = {{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};
        getAnalGateMatrix[Subscript[U|Matr|UNonNorm, __][m_]] = m;
        getAnalGateMatrix[Subscript[Ph, t__][a_]] = DiagonalMatrix[ Append[ConstantArray[1, 2^Length[{t}] - 1], Exp[I a]] ];
        getAnalGateMatrix[G[a_]] := Exp[I a] {{1,0},{0,1}};
        getAnalGateMatrix[Fac[a_]] := a {{1,0},{0,1}}; (* will not be conjugated for density matrices *)
        getAnalGateMatrix[Subscript[Rx, _][a_]] = MatrixExp[-I a/2 PauliMatrix[1]]; (* KroneckerProduct doesn't have a one-arg identity overload?? Bah *)
        getAnalGateMatrix[Subscript[Ry, _][a_]] = MatrixExp[-I a/2 PauliMatrix[2]];
        getAnalGateMatrix[Subscript[Rz, _][a_]] = MatrixExp[-I a/2 PauliMatrix[3]];
        getAnalGateMatrix[Subscript[Rx, t__][a_]] := MatrixExp[-I a/2 KroneckerProduct @@ ConstantArray[PauliMatrix[1],Length[{t}]]];
        getAnalGateMatrix[Subscript[Ry, t__][a_]] := MatrixExp[-I a/2 KroneckerProduct @@ ConstantArray[PauliMatrix[2],Length[{t}]]];
        getAnalGateMatrix[Subscript[Rz, t__][a_]] := MatrixExp[-I a/2 KroneckerProduct @@ ConstantArray[PauliMatrix[3],Length[{t}]]];
        getAnalGateMatrix[R[a_, pauli_]] := MatrixExp[-I a/2 getAnalGateMatrix @ pauli];
        getAnalGateMatrix[R[a_, paulis_Times]] := MatrixExp[-I a/2 * KroneckerProduct @@ (getAnalGateMatrix /@ List @@ paulis)]
        getAnalGateMatrix[Subscript[C, __][g_]] := getAnalGateMatrix[g]
        getAnalGateMatrix[Subscript[Id, t__]] = IdentityMatrix[2^Length[{t}]]
        
        (* extract ctrls from gate symbols *)
        getAnalGateControls[Subscript[C, c_List][___]] := c
        getAnalGateControls[Subscript[C, c__][___]] := {c}
        getAnalGateControls[_] := {}
            
        (* extract targets from gate symbols *)
        getAnalGateTargets[Subscript[U|Matr|UNonNorm, t_List][_]] := t
        getAnalGateTargets[Subscript[U|Matr|UNonNorm, t__][_]] := {t}
        getAnalGateTargets[R[_, Subscript[_, t_]]] := {t}
        getAnalGateTargets[R[_, paulis_Times]] := getAnalGateTargets /@ List @@ paulis // Flatten // Reverse
        getAnalGateTargets[Subscript[C, __][g_]] := getAnalGateTargets[g]
        getAnalGateTargets[Subscript[_, t__]] := {t}
        getAnalGateTargets[Subscript[_, t__][_]] := {t}
        
        Options[CalcCircuitMatrix] = {
            AsSuperoperator -> Automatic,
            AssertValidChannels -> True
        };
        
        (* convert a symbolic circuit channel into an analytic matrix *)
        CalcCircuitMatrix[gates_List, numQb_Integer, OptionsPattern[]] /; MemberQ[gates, Subscript[Damp|Deph|Depol|Kraus|KrausNonTP, __][__]] := 
            If[OptionValue[AsSuperoperator] =!= True && OptionValue[AsSuperoperator] =!= Automatic,
                (Message[CalcCircuitMatrix::error, "The input circuit contains decoherence channels and must be calculated as a superoperator."]; $Failed),
                With[{superops = GetCircuitSuperoperator[gates, numQb, AssertValidChannels -> OptionValue[AssertValidChannels]]},
                    If[superops === $Failed,
                    (Message[CalcCircuitMatrix::error, "Could not prepare superoperator, as per the above error."]; $Failed),
                    CalcCircuitMatrix[superops, 2*numQb]]]]
        (* convert a symbolic pure circuit into an analytic matrix *)
        CalcCircuitMatrix[gates_List, numQb_Integer, OptionsPattern[]] := 
            If[OptionValue[AsSuperoperator] === True,
                With[{superops = GetCircuitSuperoperator[gates, numQb, AssertValidChannels -> OptionValue[AssertValidChannels]]},
                    If[superops === $Failed,
                    (Message[CalcCircuitMatrix::error, "Could not prepare superoperator, as per the above error."]; $Failed),
                    CalcCircuitMatrix[superops, 2*numQb]]],
                With[{matrices = getAnalFullMatrix[
                    getAnalGateControls@#, getAnalGateTargets@#, getAnalGateMatrix@#, numQb
                    ]& /@ gates},
                    If[FreeQ[matrices, getAnalGateMatrix],
                        Dot @@ Reverse @ matrices,
                        (Message[CalcCircuitMatrix::error, "Circuit contained an unrecognised or unsupported gate: " <> 
                            ToString @ StandardForm @ First @ Cases[matrices, getAnalGateMatrix[g_] :> g, Infinity]];
                        $Failed)]]]
        CalcCircuitMatrix[gates_List, opts:OptionsPattern[]] :=
        	CalcCircuitMatrix[gates, getNumQubitsInCircuit[gates], opts]
        CalcCircuitMatrix[gate_, opts:OptionsPattern[]] :=
            CalcCircuitMatrix[{gate}, opts]
        CalcCircuitMatrix[___] := invalidArgError[CalcCircuitMatrix]
        
        GetCircuitGeneralised[gates_List] := With[
            {generalGates = Replace[gates, {
            	(* known channels are converted to Kraus maps *)
            	g:Subscript[Kraus, __][__] :> g,
            	Subscript[Damp, q_][p_] :> Subscript[Kraus, q][{
            		{{1,0},{0,Sqrt[1-p]}},
            		{{0,Sqrt[p]},{0,0}}}],
            	Subscript[Deph, q_][p_] :> Subscript[Kraus, q][{
            		Sqrt[1-p] PauliMatrix[0], 
            		Sqrt[p]   PauliMatrix[3]}],
            	Subscript[Deph, q1_,q2_][p_] :> Subscript[Kraus, q1,q2][{
            		Sqrt[1-p] IdentityMatrix[4], 
            		Sqrt[p/3] KroneckerProduct[PauliMatrix[0], PauliMatrix[3]],
            		Sqrt[p/3] KroneckerProduct[PauliMatrix[3], PauliMatrix[0]],
            		Sqrt[p/3] KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]}],
            	Subscript[Depol, q_][p_] :> Subscript[Kraus, q][{
            		Sqrt[1-p] PauliMatrix[0],
            		Sqrt[p/3] PauliMatrix[1],
            		Sqrt[p/3] PauliMatrix[2],
            		Sqrt[p/3] PauliMatrix[3]}],
            	Subscript[Depol, q1_,q2_][p_] :> Subscript[Kraus, q1,q2][ Join[
            		{Sqrt[1-p] IdentityMatrix[4]},
            		Flatten[ Table[
            			(* the PauliMatrix[0] Kraus operator is duplicated, insignificantly *)
            			If[n1===0 && n2===0, Nothing, Sqrt[p/15] KroneckerProduct[PauliMatrix[n1], PauliMatrix[n2]]],
            			{n1,{0,1,2,3}}, {n2,{0,1,2,3}}], 1]]],
            	(* global phase and fac become a fac*identity on the first qubit *)
                G[x_] :> Subscript[U, 0][ Exp[I x] IdentityMatrix[2] ],
                Fac[x_] :> Subscript[Matr, 0][ x IdentityMatrix[2] ],  (* will not be conjugated for density matrices *)
                (* Matr and UNonNorm gates remain the same *)
                g:(Subscript[Matr|UNonNorm, q__Integer|{q__Integer}][m_]) :> g,
            	(* controlled gates are turned into U of identities with bottom-right submatrix *)
                Subscript[C, c__Integer|{c__Integer}][Subscript[(mg:Matr|UNonNorm), q__Integer|{q__Integer}][m_]] :> 
                    With[{cDim=2^Length[{c}], tDim=Length@m},
                        Subscript[mg, Sequence @@ Join[{q}, {c}]][
                        Join[
                            MapThread[Join, {IdentityMatrix[cDim], ConstantArray[0,{cDim,tDim}]}],
                            MapThread[Join, {ConstantArray[0,{tDim,cDim}], m}]]]],
            	Subscript[C, c__Integer|{c__Integer}][g_] :> 
            		With[{cDim=2^Length[{c}], tDim=2^Length[getAnalGateTargets[g]]},
            			Subscript[U, Sequence @@ Join[getAnalGateTargets[g], {c}]][
            			Join[
            				MapThread[Join, {IdentityMatrix[cDim], ConstantArray[0,{cDim,tDim}]}],
            				MapThread[Join, {ConstantArray[0,{tDim,cDim}], getAnalGateMatrix@g}]]]],
            	(* all other symbols are treated like generic unitary gates *)
            	g_ :> Subscript[U, Sequence @@ getAnalGateTargets[g]][getAnalGateMatrix[g]]
            (* replace at top level *)
            }, 1]},
            If[ FreeQ[generalGates, getAnalGateMatrix],
                generalGates,
                (Message[GetCircuitGeneralised::error, "Circuit contained an unrecognised or unsupported gate: " <> 
                    ToString @ StandardForm @ First @ Cases[generalGates, getAnalGateMatrix[g_] :> g, Infinity]];
                    $Failed)
                    ]]
        GetCircuitGeneralised[op_] := GetCircuitGeneralised[{op}]
        GetCircuitGeneralised[___] := invalidArgError[GetCircuitGeneralised]
        
        getChannelAssumps[Subscript[Damp, _][x_]] := (0 <= x < 1)
        getChannelAssumps[Subscript[Deph, _][x_]] := (0 <= x < 1/2)
        getChannelAssumps[Subscript[Deph, _,_][x_]] := (0 <= x < 3/4)
        getChannelAssumps[Subscript[Depol, _][x_]] := (0 <= x < 3/4)
        getChannelAssumps[Subscript[Depol, _,_][x_]] := (0 <= x < 15/16)
        getChannelAssumps[_] := Nothing
        
        shiftInds[q__Integer|{q__Integer}, numQb_] := Sequence @@ (List[q]+numQb)
        
        Options[GetCircuitSuperoperator] = {
            AssertValidChannels -> True
        };
    
        GetCircuitSuperoperator[circ_List, numQb_Integer, OptionsPattern[]] := With[
            {superops = Flatten @ Replace[circ, {
            (* qubit-agnostic gates *)
                G[x_] :> Nothing,
                Fac[x_] :> {Fac[x]},
            (* unitaries *)
            	(* real gates (self conjugate) *)
            	Subscript[(g:H|X|Z), q__Integer|{q__Integer}] :> {Subscript[g, q], Subscript[g, shiftInds[q,numQb]]},
            	g:Subscript[P, q__Integer|{q__Integer}][v_] :> {g, Subscript[P, shiftInds[q,numQb]][v]},
            	g:Subscript[SWAP, q1_,q2_] :> {g, Subscript[SWAP, q1+numQb,q2+numQb]},
                g:Subscript[Ry, q_Integer][x_] :> {g, Subscript[Ry, q+numQb][x]},
            	(* reverse phase *)
            	Subscript[T, q_Integer] :> {Subscript[T, q], Subscript[Ph, q+numQb][-Pi/4]},
            	Subscript[S, q_Integer] :> {Subscript[S, q], Subscript[Ph, q+numQb][-Pi/2]},
            	(* reverse parameter gates (beware the mischievous Pauli Y)*)
            	Subscript[(g:Rx|Rz|Ph), q__Integer|{q__Integer}][x_] :> {Subscript[g, q][x], Subscript[g, shiftInds[q,numQb]][-x]},
                R[x_, Subscript[Y, q_Integer]] :> {R[x,Subscript[Y, q]], R[x,Subscript[Y, q+numQb]]},
            	R[x_, Subscript[p:(X|Z), q_Integer]] :> {R[x,Subscript[p, q]], R[-x,Subscript[p, q+numQb]]},
            	R[x_, Verbatim[Times][p:Subscript[_Symbol, _Integer] ..]] :> With[
                    {s = - (-1)^Count[{p}, Subscript[Y,_]]}, 
                    {R[x, Times @ p],
            		 R[s*x, Times @@ MapThread[(Subscript[#1, #2]&), {{p}[[All,1]], {p}[[All,2]] + numQb}]]}],
            	(* gates with no inbuilt conjugates *)
            	Subscript[Y, q_Integer] :> {Subscript[Y, q], Subscript[U, q+numQb][Conjugate@PauliMatrix@2]},
            	Subscript[(g:U|UNonNorm), q__Integer|{q__Integer}][m_] :> {Subscript[g, q][m], Subscript[g, shiftInds[q,numQb]][Conjugate@m]},
                Subscript[Matr, q__Integer|{q__Integer}][m_] :> {Subscript[Matr, q][m]},
            	(* controlled gates must recurse: assume inner gate resolves to two ops *)
            	Subscript[C, c__Integer|{c__Integer}][g_] :> {Subscript[C, c][g], 
            		Subscript[C, shiftInds[c,numQb]][Last @ GetCircuitSuperoperator[{g},numQb]]},
            (* channels *)
            	(* Kraus channels are turned into superoperators *)
            	Subscript[(Kraus|KrausNonTP), q__Integer|{q__Integer}][matrs_List] :> Subscript[Matr, Sequence @@ Join[{q},{shiftInds[q,numQb]}]][
            		Total[(KroneckerProduct[Conjugate[#], #]&) /@ matrs]],
            	(* other channels are first converted to Kraus, before recursing... *)
            	g:Subscript[(Damp|Depol|Deph), q__Integer|{q__Integer}][x_] :> 
            		With[{op = GetCircuitSuperoperator[GetCircuitGeneralised[g], numQb]},
                        (* and are then simplified by asserting trace-preservation *)
                        If[ OptionValue[AssertValidChannels], 
                            With[{assumps = Quiet @ Check[getChannelAssumps[g], False]},
                                (* although if CPTP is impossible, default to unsimplified, with warning *)
                                If[ assumps =!= False, 
                                    Simplify[op, assumps],
                                    Message[GetCircuitSuperoperator::error, "The channels could not be asserted as completely positive trace-preserving maps and hence were not simplified. Hide this warning with AssertValidChannels -> False, or use Quiet[]."];
                                    op]], 
                            op]],
            (* wrap unrecognised gates in dummy Head *)
            	g_ :> unrecognisedGateInSuperopCirc[g]
            (* replace at top level *)
            }, 1]},
            If[ FreeQ[superops, unrecognisedGateInSuperopCirc],
                superops,
                (Message[GetCircuitSuperoperator::error, "Circuit contained an unrecognised or unsupported gate: " <> 
                    ToString @ StandardForm @ First @ Cases[superops, unrecognisedGateInSuperopCirc[g_] :> g, Infinity]];
                    $Failed)]]
        GetCircuitSuperoperator[circ_List, opts:OptionsPattern[]] := 
            GetCircuitSuperoperator[circ, getNumQubitsInCircuit[circ], opts]
        GetCircuitSuperoperator[___] := invalidArgError[GetCircuitSuperoperator]
        
        
        
        (*
         * Below are front-end functions for 
         * modifying circuits to capture device
         * constraints and noise
         *)
         
        (* divide circ into columns, filling the left-most first *)
        GetCircuitColumns[{}] := {}
        GetCircuitColumns[circ_List] := With[{
            numQb=getNumQubitsInCircuit[circ],
            numGates=Length[circ]},
            Module[{
                gates=circ, column={}, index=1, nextIndex=Null,
                available=ConstantArray[True,numQb], compactified={}, 
                qb, i},
                
                (* continue until all gates have been grouped into a column *)
                While[index <= numGates,
                    
                    (* visit each gate from start index or until column is full (no qubits available) *)
                    For[i=index, (And[i <= numGates, Not[And @@ Not /@ available]]), i++,
                    
                        (* skip Null-marked gates (present in a previous column) *)
                        If[gates[[i]] === Null, Continue[]];
                        
                        (* extract all target and control qubits of gate (indexed from 1) *)
                        qb = 1 + Flatten @ getSymbCtrlsTargs[gates[[i]]][[{2,3}]];
                        
                        If[
                            (* if all of the gate's qubits are so far untouched in this column... *)
                            And @@ available[[qb]],
                            
                            ( (* then add the gate to the column *)
                            available[[qb]] = False;
                            AppendTo[column, gates[[i]]];
                            gates[[i]] = Null;
                            ),
                            
                            ( (* otherwise mark all the gate's qubits as unavailable (since they're blocked by this gate) *)
                            available[[qb]] = False;
                            (* and if this was the first not-in-column gate, mark for next start index *)
                            If[nextIndex === Null, nextIndex=i];
                            )
                        ]
                    ];
                    
                    (* nextIndex is unchanged if a gate occupies all qubits *)
                    If[nextIndex === Null, nextIndex=index+1];
                    
                    (* finalize the new column; empty column can arise from a trailing iteration by above edge-case *)
                    AppendTo[compactified, If[column =!= {}, column, Nothing]];
                    column = {};
                    available = ConstantArray[True,numQb];
                    index = nextIndex;
                    nextIndex = Null;
                ];
                
                (* return new gate ordering, grouped into column sub-lists *)
                compactified
            ]
        ]
            
        (* the symbolic conditions (list) under which times is monotonically increasing (and real, >0) *)
        getMotonicTimesConditions[times_List] :=
            Join[
                MapThread[Less, {Most @ times, Rest @ times} ],
                GreaterEqualThan[0] /@ times,
                (Element[#,Reals]&) /@ times]
            
        areMonotonicTimes[times_List] := 
            If[(* times must be real (else symbolic) to continue *)
                AnyTrue[(Element[#, Reals]&) /@ times, (# === False &)],
                False,
                With[{conds = getMotonicTimesConditions[times]},
                    And[
                        (* adjacent numbers increase *)
                        Not @ MemberQ[conds, False], 
                        (* symbolic numbers MAY increase (it's not impossible for them to be monotonic) *)
                        Not[False === Simplify[And @@ conds]]]]]
                        
        replaceTimeDurSymbols[expr_, spec_, timeVal_, durVal_:None] := (
            expr /. If[
                (* substitute dur value first, since it may become time dependent *)
                KeyExistsQ[spec, DurationSymbol] && durVal =!= None,
                spec[DurationSymbol] -> durVal, {}
            ] /. If[
                KeyExistsQ[spec, TimeSymbol],
                spec[TimeSymbol] -> timeVal, {}]
        )
            
        (* determines subcircuit duration, active noise and passive noise of 
         * a sub-circuit, considering the circuit variables and updating them *)
        getDurAndNoiseFromSubcirc[subcirc_, subcircTime_, spec_, forcedSubcircDur_:None] := Module[
            {activeNoises={}, passiveNoises={}, gateDurs={}, subcircDur,
             qubitActiveDurs = ConstantArray[0, spec[NumTotalQubits]], slowestGateDur},
             (* note that the final #hidden qubits of qubitActiveDurs is never non-zero *)
            
            (* iterating gates in order of appearence in subcircuit... *)
            Do[ 
                With[
                    (* determine target & control qubits *)
                    {gateQubits = Flatten @ getSymbCtrlsTargs[gate][[{2,3}]]},
                    (* attempt to match the gate against the dev spec *)
                    {gateProps = With[
                        {attemptProps = Replace[gate, spec[Gates]]},
                        Which[
                            (* if no match, throw top-level error *)
                            attemptProps === gate,
                            Throw[
                                Message[InsertCircuitNoise::error, "Encountered gate " <> ToString@StandardForm@gate <> 
                                    " which is not supported by the given device specification. Note this may be due to preceding gates," <> 
                                    " if the spec contains constraints which depend on dynamic variables. See ?GetUnsupportedGates."];
                                $Failed],
                            (* if match, but targeted qubits don't exist, throw top-level error *)
                            Not @ AllTrue[gateQubits, LessThan[spec[NumAccessibleQubits]]],
                            Throw[
                                Message[InsertCircuitNoise::error, "The gate " <> ToString@StandardForm@gate <> 
                                    " involves qubits which do not exist in the given device specification. Note that hidden qubits cannot be targeted." <> 
                                    " See ?GetUnsupportedGates."];
                                $Failed],
                            (* otherwise, gate is valid and supported by dev-spec *)
                            True,
                            attemptProps]]},
                    (* work out gate duration (time and var dependent) *)
                    {gateDur = replaceTimeDurSymbols[gateProps[GateDuration], spec, subcircTime]},    
                    (* work out active noise (time, dur and var dependent) *)
                    {gateActive = replaceTimeDurSymbols[gateProps[NoisyForm], spec, subcircTime, gateDur]},
                    (* work out var-update function (time and dur dependent) *)
                    {gateVarFunc = If[KeyExistsQ[gateProps, UpdateVariables],
                        replaceTimeDurSymbols[gateProps[UpdateVariables], spec, subcircTime, gateDur],
                        Function[None]]},
                    
                    (* collect gate info *) 
                    AppendTo[activeNoises, gateActive];
                    AppendTo[gateDurs, gateDur];
                    
                    (* update circuit variables (time, dur, var dependent, but not fixedSubcircDur dependent) *)
                    gateVarFunc[];
                    
                    (* record how long the involved qubits were activated (to later infer passive dur) *)
                    qubitActiveDurs[[ 1 + gateQubits ]] = gateDur;
                ],
                {gate,subcirc}];
                
            (* infer the whole subcircuit duration (unless forced) *)
            slowestGateDur = Max[gateDurs];
            subcircDur = If[
                forcedSubcircDur === None,
                slowestGateDur,
                forcedSubcircDur];
            (* note slowestGateDur is returned even if overriden here, for schedule-checking functions *)
                
            (* iterate all qubits (including hidden) from index 0, upward *)
            Do[
                With[
                    {qubitProps = Replace[qubit, spec[Qubits]]},
                    (* continue only if qubit matches a rule in spec[Qubits] (don't noise unspecified qubits) *)
                    If[
                        qubitProps =!= qubit,
                        With[
                            (* work out start-time and duration of qubit's passive noise (pre-determined) *)
                            {passiveTime = subcircTime + qubitActiveDurs[[ 1 + qubit ]]},
                            {passiveDur = subcircDur - qubitActiveDurs[[ 1 + qubit ]]},
                            (* work out passive noise (time, dur and var dependent *)
                            {qubitPassive = replaceTimeDurSymbols[qubitProps[PassiveNoise], spec, passiveTime, passiveDur]},
                            (* work out var-update function (time and dur dependent *)
                            {qubitVarFunc = If[KeyExistsQ[qubitProps, UpdateVariables],
                                replaceTimeDurSymbols[qubitProps[UpdateVariables], spec, passiveTime, passiveDur],
                                Function[None]]},
                                
                            (* collect qubit info *) 
                            AppendTo[passiveNoises, qubitPassive];
                            
                            (* update circuit variables (can be time, var, dur, and fixedSubcircDur dependent) *)
                            qubitVarFunc[]]]
                ],    
                {qubit, 0, spec[NumTotalQubits]-1}];
        
            (* return. Note slowestGateDur is returned even when subcircDur overrides, 
             * since that info is useful to schedule-check utilities *)
            {slowestGateDur (* subcircDur *), activeNoises, passiveNoises}
        ]
        
        getSchedAndNoiseFromSubcircs[subcircs_, spec_] := Module[
            {subcircTimes={}, subcircActives={}, subcircPassives={},
             curTime, curDur, curActive, curPassive},
                
            (* initialise circuit variables *)
            curTime = 0;
            If[KeyExistsQ[spec, InitVariables],
                spec[InitVariables][]];
            
            Do[
                (* get each subcirc's noise (updates circuit variables) *)
                {curDur, curActive, curPassive} = getDurAndNoiseFromSubcirc[sub, curTime, spec];
                    
                (* losing info here (mering separate gate infos into subcirc-wide) *)
                AppendTo[subcircTimes, curTime];
                AppendTo[subcircActives, Flatten[curActive]];
                AppendTo[subcircPassives, Flatten[curPassive]];
                
                (* keep track of the inferred schedule time *)
                curTime += curDur,
                {sub, subcircs}
            ];
            
            (* return *)    
            {subcircTimes, subcircActives, subcircPassives}
        ]
        
        getNoiseFromSched[subcircs_, subcircTimes_, spec_] := Module[
            {subcircActives={}, subcircPassives={}, subcircDurs=Differences[subcircTimes],
             dummyDur, curActive, curPassive},
             
            (* initialise circuit variables *)
            If[KeyExistsQ[spec, InitVariables],
                spec[InitVariables][]];
            
            (* if the first subcirc isn't scheduled at time=0, start with passive noise round.
             * the caller must determine this occurred (since the returned arrays are now one-item longer) *) 
            If[ First[subcircTimes] =!= 0, 
                {dummyDur, curActive, curPassive} = getDurAndNoiseFromSubcirc[{}, 0, spec, First[subcircTimes]];
                AppendTo[subcircActives, Flatten[curActive]];
                AppendTo[subcircPassives, Flatten[curPassive]];
            ];
            
            Do[
                (* get each subcirc's noise (updates circuit variables) *)
                {dummyDur, curActive, curPassive} = getDurAndNoiseFromSubcirc[subcircs[[i]], subcircTimes[[i]], spec,
                    (* force subcirc durations based on schedule, except for the final unconstrained subcirc *)
                    If[i > Length[subcircDurs], None, subcircDurs[[i]]]];
                
                AppendTo[subcircActives, Flatten[curActive]];
                AppendTo[subcircPassives, Flatten[curPassive]]; ,
                {i, Length[subcircTimes]}
            ];
            
            (* return *)    
            {subcircActives, subcircPassives}
        ]
        
        getCondsForValidSchedDurs[spec_, subcircs_, subcircTimes_] := Module[
            {forcedSubcircDurs = Differences[subcircTimes], minSubcircDur,
             dummyActive, dummyPassive, conds},
            
            (* initialise circuit variables *)
            If[KeyExistsQ[spec, InitVariables],
                spec[InitVariables][]];
            
            (* for all but the final (irrelevant) subcircuit... *)
            conds = Table[
                (* find the slowest gate (duration possibly dependent on time, 
                 * previous passive durations and vars (which are updated) *)
                {minSubcircDur, dummyActive, dummyPassive} = getDurAndNoiseFromSubcirc[ (* throws *)
                    subcircs[[i]], subcircTimes[[i]], spec, forcedSubcircDurs[[i]]];
                (* and condition that its faster than the forced subcircuit duration *)
                Function[{small, big},
                    (* if both values are numerical... *)
                    If[ NumberQ[small] && NumberQ[big],
                        (* then we need to add wiggle room for precision *)
                        Or[small <= big, Abs @ N[big-small] < 100 $MachineEpsilon],
                        (* otherwise we make a symbolic inequality *)
                        small <= big]
                ][minSubcircDur, forcedSubcircDurs[[i]]],
                {i, Length[forcedSubcircDurs]}
            ];
            
            (* by here, the validity of subcircTimes has been determined since 
             * it does not depend at all on the final subcircuit (only the durations
             * of the preceding subcircuits). However, we pass the final subcircuit 
             * to getDurAndNoiseFromSubcirc[] so that it can be checked for its general 
             * validity (gates, qubit indices, etc) and potentially throw an error. *)
             getDurAndNoiseFromSubcirc[Last@subcircs, Last@subcircTimes, spec]; (* throws *)
             
             (* return the conds *)
             conds
        ]
            
        CheckCircuitSchedule[sched:{{_, _List}..}, spec_Association] /; Not[areMonotonicTimes[sched[[All,1]]]] := (
            Message[CheckCircuitSchedule::error, "The given schedule times are not motonically increasing, nor can be for any assignment of symbols, or they are not real and positive."];
            $Failed) 
        (* this is currently naive, and assumes each sub-circuit is valid (contains simultaneous gates), and 
         * that overlapping sub-circuits is invalid. A smarter function would
         * check sub-circuits contain gates on unique qubits, and check whether 
         * overlapping sub-circuits act on unique qubits (which would be ok).
         *)
        CheckCircuitSchedule[sched:{{_, _List}..}, spec_Association] := 
            Catch[
                With[
                    {times = sched[[All,1]], subcircs = sched[[All,2]]},
                    (* list of (possibly symbolic) conditions of sufficiently long subcircuit durations *)
                    {conds = getCondsForValidSchedDurs[spec, subcircs, times]}, (* throws *)
                    (* list of (possibly symbolic) assumptions implied by monotonicity of given schedule *)
                    {mono = getMotonicTimesConditions[times]},
                    (* combine the conditions with the monotonicty constraint *)
                    With[{valid=Simplify[conds, Assumptions -> mono]},
                        Which[
                            (* if any condition is broken regardless of symbols, schedule is invalid *)
                            MemberQ[valid, False], 
                                False,
                            (* if all conditions are satisfied despite symbol assignments, the schedule is valid *)
                            AllTrue[valid, TrueQ],
                                True,
                            (* otherwise return the symbolic conditions under which the schedule is valid *)
                            True,
                                DeleteCases[valid, True]
                        ]
                    ]
                ]
            ]
        CheckCircuitSchedule[___] := invalidArgError[CheckCircuitSchedule]
        
        GetUnsupportedGates[sched:{{_, _List}..}, spec_Association] :=
            GetUnsupportedGates[ sched[[All, 2]], spec ]
        GetUnsupportedGates[cols:{_List ..}, spec_Association] :=
            GetUnsupportedGates[#, spec]& /@ cols
        GetUnsupportedGates[circ_List, spec_Association] :=
    	   Select[circ, (Not[isCompatibleGate[spec][#]]&) ]
        GetUnsupportedGates[___] := invalidArgError[GetUnsupportedGates]
            
        (* replace alias symbols (in gates & noise) with their circuit, in-place (no list nesting *)
        (* note this is overriding alias rule -> with :> which should be fine *)
        optionalReplaceAliases[False, spec_Association][in_] := in 
        (* alas, it was NOT fine; it triggered premature evaluation of the RHS of an alias rule! *)
            (*
            optionalReplaceAliases[True, spec_Association][in_] := in //. If[
                KeyExistsQ[spec, Aliases], (#1 :> Sequence @@ #2 &) @@@ spec[Aliases], {}]
            *)
        (* we will have to instead trust the user to use :> when necessary *)
        optionalReplaceAliases[True, spec_Association][in_] := in //. spec[Aliases]
        
        
        (* declaring optional args to GetCircuitSchedule *)
        Options[GetCircuitSchedule] = {
            ReplaceAliases -> False
        };
            
        (* assigns each of the given columns (unique-qubit subcircuits) a start-time *)
        GetCircuitSchedule[cols:{{__}..}, spec_Association, opts:OptionsPattern[]] := 
            Catch[
                With[
                    {times = First @ getSchedAndNoiseFromSubcircs[cols,spec]}, (* throws *)
                    Transpose[{times, cols}] // optionalReplaceAliases[OptionValue[ReplaceAliases], spec]]]
        GetCircuitSchedule[circ_List, spec_Association, opts:OptionsPattern[]] :=
            GetCircuitSchedule[GetCircuitColumns[circ], spec, opts]
        GetCircuitSchedule[___] := invalidArgError[GetCircuitSchedule]
        
        (* declaring optional args to InsertCircuitNoise *)
        Options[InsertCircuitNoise] = {
            ReplaceAliases -> False
        };

        InsertCircuitNoise[schedule:{{_, _List}..}, spec_Association, opts:OptionsPattern[]] :=
            Catch[
                (* check up-front that schedule times are monotonic *)
                If[ Not @ areMonotonicTimes @ schedule[[All,1]],
                    Throw[
                        Message[InsertCircuitNoise::error, "The scheduled subcircuit times are not monotonic (and cannot be for any substitution of any symbolic times)."];
                        $Failed]];
                        
                (* attempt to insert noise into circuit (may throw gate/circ compatibility issues) *)
                Module[
                    {times, subcircs, actives, passives},
                    {times, subcircs} = Transpose[schedule];
                    {actives, passives} = getNoiseFromSched[subcircs, times, spec]; (* throws *)
                    
                    (* pad times with initial passive noise *)
                    If[ First[times] =!= 0, PrependTo[times, 0] ];
                    (* return { {t1,subcirc1,active1,passive1}, ...} *)
                    Transpose[{times, actives, passives}] // 
                        optionalReplaceAliases[OptionValue[ReplaceAliases], spec]]]
        InsertCircuitNoise[subcircs:{{__}..}, spec_Association, opts:OptionsPattern[]] := 
            Catch[
                Transpose[getSchedAndNoiseFromSubcircs[subcircs,spec]] (* throws *) // 
                    optionalReplaceAliases[OptionValue[ReplaceAliases], spec]]   
        InsertCircuitNoise[circ_List, spec_Association, opts:OptionsPattern[]] := 
            Catch[
                Transpose[getSchedAndNoiseFromSubcircs[GetCircuitColumns[circ],spec] (*throws *) ] // 
                    optionalReplaceAliases[OptionValue[ReplaceAliases], spec]]    
        InsertCircuitNoise[___] := invalidArgError[InsertCircuitNoise]
            
        ExtractCircuit[schedule:{{_, (_List ..)}..}] :=
            Flatten @ schedule[[All,2;;]]
        ExtractCircuit[subcircs:{_List..}] :=
            Flatten @ subcircs
        ExtractCircuit[circuit_List] :=
            circuit
        ExtractCircuit[{}] := 
            {}
        ExtractCircuit[___] := invalidArgError[ExtractCircuit]
        
        formatGateParamMatrices[circ_List] :=
            circ /. m_?MatrixQ :> MatrixForm[m]
            
        ViewCircuitSchedule[sched:{{_, Repeated[_List,{1,2}]}..}, opts:OptionsPattern[]] := With[
            {isPureCirc = Length @ First @ sched == 2},
            {isActiveOnly = And[Not[isPureCirc], Flatten[ Transpose[sched][[3]] ] == {}],
             isPassiveOnly = And[Not[isPureCirc], Flatten[ Transpose[sched][[2]] ] == {}]},
            Grid[
                Which[
                    isPureCirc,
                    Join[
                        {{"time", "gates"}},
                        Function[{t,g}, {t, Row[formatGateParamMatrices[g], Spacer[0]] }] @@@ sched],
                    isActiveOnly,
                    Join[
                        {{"time", "active noise"}},
                        Function[{t,a,p}, {t, Row[formatGateParamMatrices[a], Spacer[0]] }] @@@ sched],
                    isPassiveOnly,
                    Join[
                        {{"time", "passive noise"}},
                        Function[{t,a,p}, {t, Row[formatGateParamMatrices[p], Spacer[0]] }] @@@ sched],
                    True,
                    Join[
                        {{"time", "active noise", "passive noise"}},
                        Function[{t,a,p}, {t, 
                            Row[formatGateParamMatrices[a], Spacer[0]], 
                            Row[formatGateParamMatrices[p], Spacer[0]] }] @@@ sched]
                ],    
                opts,
                Dividers -> All,
                FrameStyle -> LightGray]]
        ViewCircuitSchedule[___] := invalidArgError[ViewCircuitSchedule]
                
        (* removes suffixes $ or $123... (etc) from all symbols in expression.
           these suffixes are expected to have appeared by Mathematica's automatic 
           variable renaming in nested scoping structs (Module[ Function[]]) *)
        tidySymbolNames[exp_] :=
            exp /. s_Symbol :> RuleCondition @ Symbol @
                StringReplace[ToString[HoldForm[s]], "$"~~Repeated[NumberString,{0,1}] -> ""]
                
        (* the gates in active noise can contain symbolic qubits that won't trigger 
         * Circuit[] evaluation. This function forces Circuit[] to a list *)
        frozenCircToList[HoldForm[Circuit[gs_Times]]] := ReleaseHold[List @@@ Hold[gs]]
        frozenCircToList[HoldForm[Circuit[g_]]] := {g}
        frozenCircToList[HoldForm[gs_List]] := gs
        frozenCircToList[else_] := else
        
        (* attempt to display circuit as a column if it can be decomposed into a 
         * list (despite e.g. symbolic indices), display as is *)
        viewOperatorSeq[circ_] := With[
            {attempt = frozenCircToList[circ]},
            If[Head @ attempt === List,
                Row[attempt /. m_?MatrixQ :> MatrixForm[m], Spacer[0]],
                circ]]
                
        getHeldAssocVal[assoc_, key_] :=
            Extract[assoc, {Key[key]}, HoldForm]
        
        viewDevSpecFields[spec_, opts___] :=
            Grid[{
                {Style["Fields",Bold], SpanFromLeft},
                {"Number of accessible qubits", spec[NumAccessibleQubits]},
                {"Number of hidden qubits", spec[NumTotalQubits] - spec[NumAccessibleQubits]},
                {"Number of qubits (total)", spec[NumTotalQubits]},
                If[ KeyExistsQ[spec, TimeSymbol],
                    {"Time symbol", spec[TimeSymbol]},
                    Nothing],
                If[ KeyExistsQ[spec, DurationSymbol],
                    {"Duration symbol", spec[DurationSymbol]},
                    Nothing],
                If[KeyExistsQ[spec, InitVariables],
                    {"Variable init", spec[InitVariables]},
                    Nothing],
                {"Description", spec[DeviceDescription]}
                },
                FilterRules[{opts}, Options[Grid]],
                Dividers -> All,
                FrameStyle -> LightGray
            ] // tidySymbolNames
            
        viewDevSpecAliases[spec_, opts___] :=
            Grid[{
                {Style["Aliases",Bold], SpanFromLeft},
                {"Operator", "Definition"}
                } ~Join~ Table[
                    {
                        First[row], 
                        (* attempt to render element as spaced list *)
                        With[
                            {attemptedList = frozenCircToList[Last[row]]},
                            If[ Head[attemptedList] === List,
                                Row[attemptedList /. m_?MatrixQ :> MatrixForm[m], Spacer[0]],
                                HoldForm[attemptedList] /. m_?MatrixQ :> MatrixForm[m]
                            ]
                        ]
                    },
                    {row, List @@@ spec[Aliases]}],
                FilterRules[{opts}, Options[Grid]],
                Dividers -> All,
                FrameStyle -> LightGray
            ] // tidySymbolNames
            
        viewDevSpecActiveGates[spec_, opts___] := With[
            {showConds = Not @ FreeQ[First /@ spec[Gates], _Condition]},
            {showVars = Or @@ (KeyExistsQ[UpdateVariables] /@ Last /@ spec[Gates])},
            Grid[{
                {Style["Gates", Bold], SpanFromLeft},
                {"Gate", If[showConds,"Conditions",Nothing], "Noisy form", 
                    If[ KeyExistsQ[spec, DurationSymbol],
                        "Duration (" <> ToString@tidySymbolNames@spec[DurationSymbol] <> ")",
                        "Duration"],
                    If[showVars, "Variable update", Nothing]
                } 
                } ~Join~ Table[
                    With[
                        (*  isolate gate pattern and association *)
                        {key= First[row], props=Last[row]},
                        (* isolate gate from potential constraint *)
                        {gate=key //. c_Condition :> First[c]},
                        (* isolate potential constraints in held form *)
                        {conds=Cases[key, Verbatim[Condition][_,con_] :> HoldForm[con], {0,Infinity}, Heads -> True]},
                        {
                            gate,
                            If[showConds, Column@conds, Nothing], 
                            viewOperatorSeq @ getHeldAssocVal[props, NoisyForm], 
                            getHeldAssocVal[props, GateDuration], 
                            If[showVars, If[KeyExistsQ[props,UpdateVariables],props[UpdateVariables],""], Nothing]
                        }],
                    {row, spec[Gates]}
                ],
                FilterRules[{opts}, Options[Grid]],
                Dividers -> All,
                FrameStyle -> LightGray    
            ]  // tidySymbolNames
        ]
        
        viewDevSpecPassiveQubits[spec_, opts___] := With[
            {showVars = Or @@ (KeyExistsQ[UpdateVariables] /@ Last /@ spec[Qubits])},
            Grid[{
                {Style["Qubits", Bold], SpanFromLeft},
                {"Qubit", "Passive noise", If[showVars, "Variable update", Nothing]}
                } ~Join~ Table[
                    With[
                        {props = Replace[qubit, spec[Qubits]]},
                        {row = If[props === qubit,
                            (* show qubit with empty row for absent qubits *)
                            {qubit, "", If[showVars, "", Nothing]},
                            (* else, show its fields *)
                            {qubit,
                             viewOperatorSeq @ props[PassiveNoise] /. m_?MatrixQ :> MatrixForm[m], 
                             If[showVars, If[KeyExistsQ[props,UpdateVariables],props[UpdateVariables],""], Nothing]}]},     
                        (* insert labeled row at transition to hidden qubits *)
                        ReleaseHold @ If[qubit === spec[NumAccessibleQubits],
                            Hold[Sequence[{Style["Hidden qubits",Bold], SpanFromLeft}, row]],
                            row
                        ]
                    ],
                    {qubit, 0, spec[NumTotalQubits]-1}
                ],
                FilterRules[{opts}, Options[Grid]],
                Dividers -> All,
                FrameStyle -> LightGray 
            ] // tidySymbolNames
        ]
            
        ViewDeviceSpec[spec_Association, opts:OptionsPattern[{Grid,Column}]] := 
            Module[{view},
                Check[
                    (* get dev spec *)
                    view = Column[{
                        viewDevSpecFields[spec, opts],
                        If[KeyExistsQ[spec, Aliases] && spec[Aliases] =!= {},
                            viewDevSpecAliases[spec, opts], Nothing],
                        viewDevSpecActiveGates[spec, opts],
                        viewDevSpecPassiveQubits[spec, opts]
                        },
                        FilterRules[{opts}, Options[Column]],
                        Spacings -> {Automatic, 1}],
                    (* if errors, give warning about potential illegitimacy, and return dev spec *)
                    Echo["Note that the above errors may be illegitimate, due to premature evaluation of dynamic gate properties in the device specification. Have no fear! Device spec authors may fix this by replacing \[Rule] with \[RuleDelayed] in NoisyForm and GateDuration of Gates which feature variables.", "ViewDeviceSpec: "];
                    view
                ]
            ]    
        ViewDeviceSpec[___] := invalidArgError[ViewDeviceSpec]
        
        getDeviceSpecIssueString[spec_] := Catch[
            
            (* check top-level required keys *)
            Do[
                If[ Not @ KeyExistsQ[spec, key], 
                    Throw["Specification is missing the required key: " <> SymbolName[key] <> "."]],
                {key, {DeviceDescription, NumAccessibleQubits, NumTotalQubits, 
                       Gates, Qubits}}];
            
            (* check number of qubits *)
            Do[ 
                If[ Not @ MatchQ[spec[key], n_Integer /; n > 0 ], 
                    Throw["NumAccessibleQubits and NumTotalQubits must be positive integers."]],
                {key, {NumAccessibleQubits, NumTotalQubits}}];

            If[ spec[NumAccessibleQubits] > spec[NumTotalQubits],
                Throw["NumAccessibleQubits cannot exceed NumTotalQubits."]];
            
            (* check symbols are indeed symbols *)
            Do[
                If[ KeyExistsQ[spec, key],
                    If[ Not @ MatchQ[spec[key], _Symbol], 
                        Throw["TimeSymbol and DurationSymbol must be symbols."] ]],
                {key, {TimeSymbol, DurationSymbol}}];
                
            (* check  alias is a list of delayed rules to circuits *)
            If[ KeyExistsQ[spec, Aliases], 
                If[ Not @ MatchQ[ spec[Aliases], { (_ :>  (_Circuit | _List)) ... } ],
                    Throw["Aliases must be a list of DelayedRule, each pointing to a Circuit (or a list of operators)."]]]
                
            (* check aliases do not contain symbols *)
            If[ KeyExistsQ[spec, Aliases], 
                Do[
                    If[ KeyExistsQ[spec, key],
                        If[ Not @ FreeQ[spec[Aliases], spec[key]], 
                            Throw["Aliases (definitions or operators) must not feature TimeSymbol nor DurationSymbol; they can instead be passed as arguments to the alias operator."]]],
                    {key, {TimeSymbol, DurationSymbol}}]];
                    
            (* check alias LHS don't include conditions *)
            If[ KeyExistsQ[spec, Aliases], 
                If[ Not @ FreeQ[First /@ spec[Aliases], Condition],
                    Throw["Aliases must not include Condition in their operators (the left-hand side of RuleDelayed)."]]];
                
            (* check init-var is zero-arg function *)
            If[ KeyExistsQ[spec, InitVariables],
                If[ Not[
                    MatchQ[ spec[InitVariables], _Function ] &&
                    Quiet @ Check[ spec[InitVariables][]; True, False ] ], (* duck typed *)
                    Throw["InitVariables must be a zero-argument Function (or excluded entirely), which initialises any variables needing later modification in UpdateVariables."]]];
            
            (* check Gates and Qubits are list of RuleDelayed, to an association *)
            Do[
                (* If[ (Not @ MatchQ[spec[key], _List]) || (Not @ AllTrue[spec[key], MatchQ[RuleDelayed[_, _Association]]]), *)
                If[ Not @ MatchQ[spec[key], { (_ :>  _Association) ... }],
                    Throw["Gates and Qubits must each be a list of RuleDelayed, each pointing to an Association."]],
                {key, {Gates,Qubits}}];
                
            (* check every Gates association has required keys *)
            Do[
                If[ Not @ KeyExistsQ[assoc, key],
                    Throw["An Association in Gates is missing required key " <> SymbolName[key] <> "."]],
                {assoc, Last /@ spec[Gates]},
                {key, {NoisyForm, GateDuration}}];
                
            (* check that Gates patterns do not refer to symbols *)
            Do[
                If[ KeyExistsQ[spec, key],
                    If[ Not @ FreeQ[pattern, spec[key]],
                        Throw["The operator patterns in Gates (left-hand side of the rules) must not include the TimeSymbol or the DurationSymbol (though the right-hand side may)."]]],
                {pattern, First /@ spec[Gates]},
                {key, {TimeSymbol, DurationSymbol}}];
                
            (* check every Gates' GateDuration doesn't contain the duration symbol (self-reference) *)
            If[ KeyExistsQ[spec, DurationSymbol],  
                Do[
                    If[ Not @ FreeQ[dur, spec[DurationSymbol]],
                        Throw["A GateDuration cannot refer to the DurationSymbol, since the DurationSymbol is substituted the value of the former."]],
                    {dur, spec[Gates][[All, 2]][GateDuration] // Through}]];
            
            (* check every Qubit assoc contains required keys *)
            Do[
                If[ Not @ KeyExistsQ[assoc, PassiveNoise],
                    Throw["An association in Qubits is missing required key PassiveNoise."]],
                {assoc, Last /@ spec[Qubits]}]
                
            (* check every passive noise is a list (or Circuit, not yet evaluating) *)
            Do[
                If[ Not @ MatchQ[passive, _List|_Circuit],
                    Throw["Each PassiveNoise must be a Circuit[] or list of operators."]],
                {passive, spec[Qubits][[All, 2]][PassiveNoise] // Through}];
                
            (* check every update-vars (in Gates and Qubits assoc) is a zero-arg function *)
            Do[
                If[ KeyExistsQ[assoc, UpdateVariables],
                    If[ Not[
                        MatchQ[ assoc[UpdateVariables], _Function ] &&
                        (* calling UpdateVariables[] may generate other errors, but we care only about zero-args *)
                        Quiet @ Check[ assoc[UpdateVariables][]; True, False, Function::fpct ] ], (* duck typed *)
                        Throw["Each UpdateVariables must be a zero-argument Function (or excluded entirely)."]]],
                {assoc, Join[spec[Gates][[All,2]], spec[Qubits][[All,2]]] }];
                
            (* no detected issues *)
            None
        ]
        
        CheckDeviceSpec[spec_Association] := With[
            {issue = Quiet @ getDeviceSpecIssueString[spec]},
            If[ issue === None, 
                True,
                Message[CheckDeviceSpec::error, issue]; False]]    
        CheckDeviceSpec[___] := (
            Message[CheckDeviceSpec::error, "Argument must be a single Association."];
            $Failed)
            
        
        
        (*
         * Below are front-end functions 
         * for modifying circuits
         *)    
        
        getInverseGate[g:Subscript[H|Id|SWAP|X|Y|Z, __]] := g
        getInverseGate[(g:Subscript[Rx|Ry|Rz|Ph, __])[x_]] := g[-x]
        getInverseGate[R[x_, s_]] := R[-x, s]
        getInverseGate[Subscript[T, q_]] := Subscript[Ph, q][-Pi/4]
        getInverseGate[Subscript[S, q_]] := Subscript[Ph, q][-Pi/2]
        getInverseGate[G[x_]] := G[-x]
        getInverseGate[Fac[x_]] := Fac[1/x]
        getInverseGate[Subscript[(g:U|UNonNorm), q__][m_?MatrixQ]] := Subscript[g, q][ConjugateTranspose[m]]
        getInverseGate[Subscript[Matr, q__][m_?MatrixQ]] := Subscript[Matr, q][Inverse[m]]
        getInverseGate[g:Subscript[C, c__][h_]] := With[
            {hInv = getInverseGate[h]},
            If[Head @ hInv =!= $Failed, 
                Subscript[C, c][hInv], 
                $Failed[g]]]
        getInverseGate[g_] := $Failed[g]

        GetCircuitInverse[circ_List] := With[
            {invs = getInverseGate /@ Reverse[circ]},
            {bad = FirstCase[invs, $Failed[g_] :> g]},
            If[bad === Missing["NotFound"],
                invs,
                (Message[GetCircuitInverse::error, "Could not determine the inverse of gate " <> ToString@TraditionalForm@bad <> "."];
                $Failed)]]
        GetCircuitInverse[___] := invalidArgError[GetCircuitInverse]
        
        tidyInds[q__] := Sequence @@ Sort@DeleteDuplicates@List@q
        tidyMatrixGate[Subscript[g:(U|Matr|UNonNorm), q_Integer][m_]] := Subscript[g, q][Simplify @ m]
        tidyMatrixGate[Subscript[g:(U|Matr|UNonNorm), q__Integer][m_]] /; OrderedQ[{q}] := Subscript[g, q][Simplify @ m]
        tidyMatrixGate[Subscript[g:(U|Matr|UNonNorm), q__Integer][m_]] := 
        	With[{order=Ordering[{q}]},
        		Do[
        			If[order[[i]] =!= i, Block[{tmp}, With[
        				{q1={q}[[i]], q2={q}[[order[[i]]]]}, 
        				{s=CalcCircuitMatrix[{Subscript[SWAP, i-1,order[[i]]-1]}, Length[{q}]]},
        				Return @ tidyMatrixGate @ Subscript[g, Sequence@@((({q} /. q1->tmp) /. q2->q1) /. tmp->q2)][s . m . s]]]],
        			{i, Length[{q}]}]]

        SimplifyCircuit[circ_List] := With[{
        	(*
        	 * establish preconditions
        	 *)
        	initCols = GetCircuitColumns[circ] //. {
        		(* convert Ph controls into targets *)
        		Subscript[C, c__Integer|{c__Integer}][Subscript[Ph, t__Integer|{t__Integer}][x__]] :> Subscript[Ph, tidyInds[c,t]][x],
        		(* convert S and T gates into Ph *)
        		Subscript[(g:S|T), t_Integer ]:> Subscript[Ph, t][Pi / (g/.{S->2,T->4})], 
        		Subscript[C, c__Integer|{c__Integer}][Subscript[(g:S|T), t_Integer]] :> Subscript[Ph, tidyInds[c,t]][Pi / (g/.{S->2,T->4})],
        		(* sort qubits of general unitaries by SWAPs upon matrix *)
        		g:Subscript[U|Matr|UNonNorm, q__Integer|{q__Integer}][m_] :> tidyMatrixGate[g],
        		Subscript[C, c__][g:Subscript[U|Matr|UNonNorm, q__Integer|{q__Integer}][m_]] :> Subscript[C, tidyInds@c][tidyMatrixGate@g],
        		(* sort controls of any gate *)
        		Subscript[C, c__Integer|{c__Integer}][g_] :> Subscript[C, tidyInds@c][g],
        		(* sort targets of target-order-agnostic gates *)
        		Subscript[(g:(H|X|Y|Z|Id|SWAP|Ph|M||T|S)), t__Integer|{t__Integer}] :> Subscript[g, tidyInds@t],
        		Subscript[(g:(Rx|Ry|Rz|Damp|Deph|Depol)), t__Integer|{t__Integer}][x__] :> Subscript[g, tidyInds@t][x],
        		(* unpack all qubit lists *)
        		Subscript[s_, {t__}] :> Subscript[s, t],
        		Subscript[s_, {t__}][x_] :> Subscript[s, t][x],
        		Subscript[C, c__][Subscript[s_, {t__}]] :> Subscript[C, c][Subscript[s, t]],
        		Subscript[C, c__][Subscript[s_, {t__}][x_]] :> Subscript[C, c][Subscript[s, t][x]]
        	}},
        	(* above establishes preconditions:
        		- gates within a column target unique qubits
        		- qubit lists of order-agnostic gates are ordered and duplicate-free
        		- qubit lists are flat (not contained in List)
        		- phase gates have no control qubits
        		- there are no T or S gates
        		- R[x, pauli-tensor] have fixed-order tensors (automatic by Times)
        		- the first global phase G[] will appear in the first column
        	*)
        	Module[{simpCols},
        		(* 
        		* repeatedly simplify circuit until static 
        		*)
        		simpCols = FixedPoint[ Function[{prevCols}, Module[{cols},
        			cols = prevCols;
        			(* 
        			 * simplify contiguous columns
        			 *)
        			cols = SequenceReplace[cols, Join[
        				(* remove adjacent idempotent operations *)
        				Join @@ Table[
        					{ {a___, wrap@Subscript[gate, q__], b___}, {c___,  wrap@Subscript[gate, q__], d___} } :> Sequence[{a,b},{c,d}],
        					{gate, {H,X,Y,Z,Id,SWAP}},
        					{wrap, {Identity, Subscript[C, ctrls__]}}],
        				(* combine arguments of adjacent parameterized gates *)
        				Join @@ Table[ 
        					(* awkward With[] use to force immediate eval of 'gate' in DelayedRule *)
        					With[{gate=gateSymb}, {
        					{ {a___, Subscript[gate, q__][x_], b___}, {c___, Subscript[gate, q__][y_], d___} } :> Sequence[{a,Subscript[gate, q][x+y//Simplify],b},{c,d}],
        					{ {a___, Subscript[C, ctrl__][Subscript[gate, q__][x_]], b___}, {c___, Subscript[C, ctrl__][Subscript[gate, q__][y_]], d___} } :> Sequence[{a,Subscript[C, ctrl][Subscript[gate, q][x+y//Simplify]],b},{c,d}]
        					}],
        					{gateSymb, {Ph,Rx,Ry,Rz}}],
        				{
        					{ {a___, R[x_,op_], b___}, {c___, R[y_,op_], d___} } :> Sequence[{a,R[x+y//Simplify,op],b},{c,d}],
        					{ {a___, Subscript[C, ctrl__]@R[x_,op_], b___}, {c___, Subscript[C, ctrl__]@R[y_,op_], d___} } :> Sequence[{a,Subscript[C, ctrl]@R[x+y//Simplify,op],b},{c,d}]
        				},
        				(* multiply matrices of adjacent unitaries and Matr *)
                        (* we do not presently merge U unto neighbouring Matr *)
        				{
        					{ {a___, Subscript[g:(U|Matr|UNonNorm), q__][m1_], b___}, {c___, Subscript[g:(U|Matr|UNonNorm), q__][m2_], d___} } :> Sequence[{a,Subscript[g, q][m1 . m2//Simplify],b},{c,d}],
        					{ {a___, Subscript[C, ctrl__]@Subscript[g:(U|Matr|UNonNorm), q__][m1_], b___}, {c___, Subscript[C, ctrl__]@Subscript[g:(U|Matr|UNonNorm), q__][m2_], d___} } :> Sequence[{a,Subscript[C, ctrl]@Subscript[g, q][m1 . m2//Simplify],b},{c,d}]
        				},
        				(* merge all global phases *)
        				{
        					{ {a___, G[x_], b___, G[y_], c___} } :> {G[x+y//Simplify], a, b, c},
        					{ {a___, G[x_], b___}, infix___, {c___, G[y_], d___} } :> Sequence[{G[x+y//Simplify],a,b}, infix, {c,d}]
        				},
                        (* merge all factors *)
                        {
        					{ {a___, Fac[x_], b___, Fac[y_], c___} } :> {Fac[x y //Simplify], a, b, c},
        					{ {a___, Fac[x_], b___}, infix___, {c___, Fac[y_], d___} } :> Sequence[{Fac[x y //Simplify],a,b}, infix, {c,d}]
        				},
                        (* merge factors and global phases *)
                        {
                            { {a___, Fac[x_], b___, G[y_], c___} } :> {Fac[x Exp[y I] //Simplify], a, b, c},
                            { {a___, G[y_], b___, Fac[x_], c___} } :> {Fac[x Exp[y I] //Simplify], a, b, c},
        					{ {a___, Fac[x_], b___}, infix___, {c___, G[y_], d___} } :> Sequence[{Fac[x Exp[y I] //Simplify],a,b}, infix, {c,d}],
                            { {a___, G[y_], b___}, infix___, {c___, Fac[x_], d___} } :> Sequence[{Fac[x Exp[y I] //Simplify],a,b}, infix, {c,d}]
                        },
        				(* merge adjacent Pauli operators *)
        				{
        					{ {a___, Subscript[X, q__], b___}, {c___, Subscript[Y, q__], d___} } :> Sequence[{a,G[3 Pi/2],Subscript[Z, q],b},{c,d}],
        					{ {a___, Subscript[Y, q__], b___}, {c___, Subscript[Z, q__], d___} } :> Sequence[{a,G[3 Pi/2],Subscript[X, q],b},{c,d}],
        					{ {a___, Subscript[Z, q__], b___}, {c___, Subscript[X, q__], d___} } :> Sequence[{a,G[3 Pi/2],Subscript[Y, q],b},{c,d}],
        					{ {a___, Subscript[Y, q__], b___}, {c___, Subscript[X, q__], d___} } :> Sequence[{a,G[Pi/2],Subscript[Z, q],b},{c,d}],
        					{ {a___, Subscript[Z, q__], b___}, {c___, Subscript[Y, q__], d___} } :> Sequence[{a,G[Pi/2],Subscript[X, q],b},{c,d}],
        					{ {a___, Subscript[X, q__], b___}, {c___, Subscript[Z, q__], d___} } :> Sequence[{a,G[Pi/2],Subscript[Y, q],b},{c,d}]
        				},
        				(* merge adjacent rotations with paulis *)
        				Join @@ Table[ With[{rot=First@ops, pauli=Last@ops}, {
        					{ {a___, Subscript[pauli, q__], b___}, {c___, Subscript[rot, q__][x_], d___} } :> Sequence[{a,G[Pi/2],Subscript[rot, q][x+Pi],b},{c,d}],
        					{ {a___, Subscript[rot, q__][x_], b___}, {c___, Subscript[pauli, q__], d___} } :> Sequence[{a,G[Pi/2],Subscript[rot, q][x+Pi],b},{c,d}]
        				}], 
        					{ops, {{Rx,X},{Ry,Y},{Rz,Z}}}]
        				
        				(* TODO: should I convert Z to Ph too in order to compound with Ph?? *)
        				(* and should controlled Z become Ph too?? *)
        				
        				(* TODO: turn Rx[2\[Pi] + eh] = G[\[Pi]] Rx[eh] ??? *)
        			]];
        			(* 
        			 * simplify single gates (at any level, even within controls)
        			 *)
        			cols = cols //. {
        				(* remove empty columns *)
        				{} -> Nothing, 
        				{{}..} -> Nothing,
        				(* remove controlled gates with insides already removed *)
        				Subscript[C, __][Nothing] -> Nothing,
        				(* remove zero-parameter gates *)
        				Subscript[(Ph|Rx|Ry|Rz|Damp|Deph|Depol), __][0] -> Nothing,
        				R[0,_] -> Nothing,
        				G[0] -> Nothing,
                        Fac[1|1.] -> Nothing,
                        Fac[x_ /; (Abs[x] === 1)] -> G[ ArcTan[Re@x, Im@x] ],
        				(* remove identity matrices (qubits are sorted) *)
        				Subscript[U|Matr|UNonNorm, q__][m_] /; m === IdentityMatrix[2^Length[{q}]] -> Nothing,
        				(* simplify known parameters to within their periods *)
        				Subscript[Ph, q__][x_?NumericQ] /; Not[0 <= x < 2 Pi] :> Subscript[Ph, q]@Mod[x, 2 Pi],
        				(g:(Subscript[(Rx|Ry|Rz), q__]))[x_?NumericQ] /; Not[0 <= x < 4 Pi] :> g@Mod[x, 4 Pi],
        				R[x_?NumericQ, op_] /; Not[0 <= x < 4 Pi] :> R[Mod[x, 4 Pi], op],
        				(* convert single-target R to Rx, Ry or Rz *)
        				R[x_, op:Subscript[(X|Y|Z), q_]] :> (op[x] /. {X->Rx,Y->Ry,Z->Rz})
        			};
        			(* 
        			 * simplify single gates (top-level only, cannot occur within controls) 
        			 *)
        			cols = Replace[cols, {
        				(* transform param gates with Pi params *)
        				(g:(Subscript[(Rx|Ry|Rz), q__]))[Pi] :> Sequence[ G[3 Pi/2], g/.{Rx->X,Ry->Y,Rz->Z} ],
        				R[Pi, op_Times] :> Sequence[G[3 Pi/2], Sequence@@op],
        				Subscript[Ph, q_][Pi] :> Subscript[Z, q],
        				Subscript[Ph, q__][Pi] :> Subscript[C, Sequence@@Rest[{q}]][Subscript[Z, First@{q}]],
        				(* transform rotations gates with 2 Pi params *)
        				(g:(Subscript[(Rx|Ry|Rz), q__]))[2 Pi] :> G[Pi],
        				R[2 Pi, _Times] :> G[Pi]
        				(* forces top-level, inside each subcircuit *)
        				}, {2}];
        			(* 
        			 * re-update circuit columns 
        			 *)
        			If[cols === Nothing, cols={}];
        			cols = GetCircuitColumns @ ExtractCircuit @ cols;
        			cols
        		]], initCols];
        		(*
        		 * post-process the simplified columns
        		 *)
        		ExtractCircuit @ simpCols /. {
        			Subscript[Ph, q_][Pi/2] :> Subscript[S, q],
        			Subscript[Ph, q_][Pi/4] :> Subscript[T, q]
        			
        			(* TODO: controls?? Z gates?? *)
        		}]]
        
        SimplifyCircuit[___] := invalidArgError[SimplifyCircuit]
        
        
        
        (*
         * Below are front-end functions 
         * for generating circuits for 
         * GetKnownCircuit[]
         *)
         
        GetKnownCircuit["QFT", qubits_List] := Flatten @ {
            Table[ { 
                Subscript[H, qubits[[n]]], 
                Table[Subscript[Ph, qubits[[n]],qubits[[n-m]]][Pi/2^m], {m,1,n-1}]},
                {n, Length[qubits], 1, -1}],
            Table[
                Subscript[SWAP, qubits[[q]], qubits[[-q]]], 
                {q, 1, Floor[Length[qubits]/2]}] }
        GetKnownCircuit["QFT", numQubits_Integer] :=
            GetKnownCircuit["QFT", Range[0,numQubits-1]]
            
        separateCoeffAndPauliTensor[pauli_Subscript] := {1, pauli}
        separateCoeffAndPauliTensor[prod_Times] := {
        	Times@@Cases[prod, c:Except[Subscript[(X|Y|Z|Id), _Integer]]:>c],
        	Times@@Cases[prod, p:Subscript[(X|Y|Z|Id), _Integer]:>p] }
        separateTermsOfPauliHamil[hamil_Plus] := 
        	separateCoeffAndPauliTensor /@ (List @@ hamil)
        separateTermsOfPauliHamil[term_] := 
        	{separateCoeffAndPauliTensor[term]}
        getSymmetrizedTerms[terms_List, fac_, 1] := 
        	MapAt[fac # &, terms, {All, 1}]
        getSymmetrizedTerms[terms_List, fac_, 2] := With[
        	{s1 = getSymmetrizedTerms[terms, fac/2, 1]}, 
        	Join[Most[s1], {{2 s1[[-1,1]], s1[[-1,2]]}}, Reverse[Most[s1]]]]
        getSymmetrizedTerms[terms_List, fac_, n_?EvenQ] := 
        	Block[{x, p=1/(4-4^(1/(n-1)))}, With[
        		{s = getSymmetrizedTerms[terms, x, n-2]}, 
        		{r = s /. x -> fac p},
        		Join[r, r, s /. x -> (1-4p)fac, r, r]]]
        getTrotterTerms[terms_List, order_, reps_, time_] :=
        	With[{s=getSymmetrizedTerms[terms, time/reps, order]},
        		Join @@ ConstantArray[s, reps]]
                
        GetKnownCircuit["Trotter", hamil_, order_Integer, reps_Integer, time_] /; (
        	order>=1 && (order===1 || EvenQ@order) && reps>=1) := 
        	With[
        		{terms = separateTermsOfPauliHamil @ hamil},
        		{gates = (R[2 #1, #2]&) @@@ getTrotterTerms[terms, order, reps, time]},
        		gates /. R[_, Subscript[Id, _Integer]] :> Nothing]
                
        GetKnownCircuit["HardwareEfficientAnsatz", reps_Integer, param_Symbol, qubits_List] := 
        	Module[{i=1, ent},
                ent = Subscript[C, #[[1]]][Subscript[Z, #[[2]]]]& /@ Partition[qubits,2,1,{1,1}];
            	Flatten[{
            		Table[{
            			Table[Subscript[g, q][param[i++]], {q,qubits}, {g,{Ry,Rz}}],
            			ent[[1 ;; ;; 2 ]],
                        ent[[2 ;; ;; 2 ]]
            		}, reps],
            		Table[Subscript[g, q][param[i++]], {q,qubits}, {g,{Ry,Rz}}]}]]
        GetKnownCircuit["HardwareEfficientAnsatz", reps_Integer, param_Symbol, numQubits_Integer] :=
        	GetKnownCircuit["HardwareEfficientAnsatz", reps, param, Range[0,numQubits-1]]
            
        GetKnownCircuit["TrotterAnsatz", hamil_, order_Integer, reps_Integer, param_Symbol] := 
            Module[{i=1},
                GetKnownCircuit["Trotter", hamil, order, reps, 1] /. {
                    R[x_, p_] :> R[param[i++], p],
                    (g:Subscript[Rx|Ry|Rz, _])[x_] :> g[param[i++]]}]
                    
        GetKnownCircuit["LowDepthAnsatz", reps_Integer, paramSymbol_Symbol, qubits_List] := 
            Module[{i=1, pairs=Most@Partition[qubits,2,1,{1,1}]},
                Flatten @ Join[
                    Table[Subscript[Rz, q][paramSymbol[i++]], {q,qubits}],
                    Table[{
                            R[ paramSymbol[i++], Subscript[X, #1] Subscript[Y, #2] ],
                            R[ paramSymbol[i++], Subscript[Y, #1] Subscript[X, #2] ],
                            R[ paramSymbol[i++], Subscript[Y, #1] Subscript[Y, #2] ],
                            R[ paramSymbol[i++], Subscript[X, #1] Subscript[X, #2] ]
                            }& @@@
                            Join[ pairs[[1;;;;2]], pairs[[2;;;;2]] ],
                        reps]]]
        GetKnownCircuit["LowDepthAnsatz", reps_Integer, paramSymbol_Symbol, numQubits_Integer] :=
            GetKnownCircuit["LowDepthAnsatz", reps, paramSymbol, Range[0,numQubits-1]]
        GetKnownCircuit[___] := invalidArgError[GetKnownCircuit]

    End[ ]
                                       
EndPackage[]

Needs["QuEST`Option`"]

Needs["QuEST`Gate`"]

Needs["QuEST`DeviceSpec`"]

Needs["QuEST`Deprecated`"]
