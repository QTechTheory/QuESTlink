
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
    
    CalcQuregDerivs::usage = "CalcQuregDerivs[circuit, initQureg, varVals, derivQuregs] sets the given list of (deriv)quregs to be the result of applying derivatives of the parameterised circuit to the initial state. The derivQuregs are ordered by the varVals, which should be in the format {param -> value}, where param is featured in Rx, Ry, Rz, R or U (and controlled) of the given circuit ONCE (multiple times within a U matrix is allowed). The initQureg is unchanged. Note Rx[theta] is allowed, but Rx[f(theta)] is not. Furthermore U matrices must contain at most one parameter."
    CalcQuregDerivs::error = "`1`"
    
    CalcInnerProducts::usage = "CalcInnerProducts[quregIds] returns a Hermitian matrix with i-th j-th element CalcInnerProduct[quregIds[i], quregIds[j]].
CalcInnerProducts[braId, ketIds] returns a complex vector with i-th element CalcInnerProduct[braId, ketIds[i]]."
    CalcInnerProducts::error = "`1`"

    CalcDensityInnerProducts::usage = "CalcDensityInnerProducts[quregIds] returns a real, symmetric matrix with i-th j-th element CalcDensityInnerProduct[quregIds[i], quregIds[j]].
CalcDensityInnerProducts[rhoId, omegaIds] returns a real vector with i-th element CalcDensityInnerProduct[rhoId, omegaIds[i]]."
    CalcDensityInnerProducts::error = "`1`"
    
    Circuit::usage = "Circuit[gates] converts a product of gates into a left-to-right circuit, preserving order."
    Circuit::error = "`1`"
    
    Operator::usage = "Operator[gates] converts a product of gates into a right-to-left circuit."
    Operator::error = "`1`"
    
    CalcExpecPauliProd::usage = "CalcExpecPauliProd[qureg, paulis, workspace] evaluates the expected value of a product of Paulis. workspace must be a qureg of equal dimensions to qureg."
    CalcExpecPauliProd::error = "`1`"

    CalcExpecPauliSum::usage = "CalcExpecPauliSum[qureg, pauliSum, workspace] evaluates the expected value of a weighted sum of Pauli products, of a normalised qureg. workspace must be a qureg of equal dimensions to qureg. qureg is unchanged, and workspace is modified."
    CalcExpecPauliSum::error = "`1`"

    ApplyPauliSum::usage = "ApplyPauliSum[inQureg, pauliSum, outQureg] modifies outQureg to be the result of applying the weighted sum of Paulis to inQureg."
    ApplyPauliSum::error = "`1`"
    
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
    
    CalcPauliSumMatrix::usage = "CalcPauliSumMatrix[pauliSum] returns the numerical matrix of the given real-weighted sum of Pauli operators. The number of qubits is assumed to be the largest Pauli target. This accepts only sums of Pauli products with unique qubits and floating-point coefficients, and is computed numerically."
    CalcPauliSumMatrix::error = "`1`"
    
    CalcPauliExpressionMatrix::usage = "CalcPauliExpressionMatrix[expr] returns the analytic matrix given by the symbolic expression of Pauli operators, X, Y, Z, Id. The number of qubits is assumed to be the largest Pauli target. Accepts the same inputs as SimplfyPaulis[], and is computed symbolically"
    CalcPauliExpressionMatrix::error = "`1`"

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
    
    GetPauliSumFromCoeffs::usage = "GetPauliSumFromCoeffs[addr] opens or downloads the file at addr (a string, of a file location or URL), and interprets it as a list of coefficients and Pauli codes, converting this to a symbolic weighted sum of Pauli products. Each line of the file is a separate term (a Pauli product), with format {coeff code1 code2 ... codeN} (exclude braces) where the codes are in {0,1,2,3} (indicating a I, X, Y, Z term in the product respectively), for an N-qubit operator. Each line must have N+1 terms (including the real decimal coefficient at the beginning)."
    GetPauliSumFromCoeffs::error = "`1`"
    
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

    CalcCircuitMatrix::usage = "CalcCircuitMatrix[circuit] returns an analytic expression for the given unitary circuit, which may contain undefined symbols. The number of qubits is inferred from the circuit indices (0 to maximum specified).
CalcCircuitMatrix[circuit, numQubits] gives CalcCircuitMatrix a clue about the number of present qubits."
    CalcCircuitMatrix::error = "`1`"
    
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
    
    ExtractCircuit::usage = "ExtractCircuit[] returns the ultimate circuit from the outputs of InsertCircuitNoise[], GetCircuitSchedule[] and GetCircuitSchedule[]."
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
GetKnownCircuit[\"Trotter\", hamil, order, reps, time]"
    GetKnownCircuit::error = "`1`"
    
    
    (*
     * optional arguments to public functions
     *)
     
    BeginPackage["`Option`"]

    WithBackup::usage = "Optional argument to ApplyCircuit, indicating whether to create a backup during circuit evaluation to restore the input state in case of a circuit error. This incurs additional memory (default True). If the circuit contains no error, this option has no effect besides wasting memory."
    
    ShowProgress::usage = "Optional argument to ApplyCircuit, indicating whether to show a progress bar during circuit evaluation (default False). This slows evaluation slightly."
    
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
    
    EndPackage[]
    
    
    
    (* 
     * gate symbols    
     *)
     
    BeginPackage["`Gate`"]

    H::usage = "H is the Hadamard gate."
    
    X::usage = "X is the Pauli X gate, a.k.a NOT or bit-flip gate."
    
    Y::usage = "Y is the Pauli Y gate."
    
    Z::usage = "Z is the Pauli Z gate."
    
    Rx::usage = "Rx[\[Theta]] is a rotation of \[Theta] around the x-axis of the Bloch sphere, Exp[-\[ImaginaryI] \[Theta]/2 X \[CircleTimes] X \[CircleTimes]...]."        
    
    Ry::usage = "Ry[\[Theta]] is a rotation of \[Theta] around the y-axis of the Bloch sphere, Exp[-\[ImaginaryI] \[Theta]/2 Y \[CircleTimes] Y \[CircleTimes]...]." 
    
    Rz::usage = "Rz[\[Theta]] is a rotation of \[Theta] around the z-axis of the Bloch sphere, Exp[-\[ImaginaryI] \[Theta]/2 Z \[CircleTimes] Z \[CircleTimes]...]." 
    
    R::usage = "R[\[Theta], paulis] is the unitary Exp[-\[ImaginaryI] \[Theta]/2 \[CircleTimes] paulis]."   
    
    S::usage = "S is the S gate, a.k.a. PI/2 gate."
    
    T::usage = "T is the T gate, a.k.a PI/4 gate."
    
    U::usage = "U[matrix] is a general 1 or 2 qubit unitary gate, enacting the given 2x2 or 4x4 matrix."
    
    Deph::usage = "Deph[prob] is a 1 or 2 qubit dephasing with probability prob of error."
    
    Depol::usage = "Depol[prob] is a 1 or 2 qubit depolarising with probability prob of error."
    
    Damp::usage = "Damp[prob] is 1 qubit amplitude damping with the given decay probability."
    
    SWAP::usage = "SWAP is a 2 qubit gate which swaps the state of two qubits."
    
    M::usage = "M is a destructive measurement gate which measures the indicated qubits in the Z basis. Targeting multiple qubits is the same as applying M to each in-turn, though their outcomes will be grouped in the output of ApplyCircit[]."
    
    P::usage = "P[val] is a (normalised) projector onto {0,1} (i.e. a forced measurement) such that the target qubits represent integer val in binary (right most target takes the least significant digit in val).
P[outcome1, outcome2, ...] is a (normalised) projector onto the given {0,1} outcomes. The left most qubit is set to the left most outcome.
The probability of the forced measurement outcome (if hypothetically not forced) is included in the output of ApplyCircuit[]."
    
    Kraus::usage = "Kraus[ops] applies a one or two-qubit Kraus map (given as a list of Kraus operators) to a density matrix."
    
    KrausNonTP::usage = "Kraus[ops] applies a one or two-qubit non-trace-preserving Kraus map (given as a list of matrix operators) to a density matrix."
    
    G::usage = "G[\[Theta]] applies a global phase rotation of phi, by premultiplying Exp[\[ImaginaryI] \[Theta]]."
    
    Id::usage = "Id is an identity gate which effects no change, but can be used for forcing gate alignment in DrawCircuit, or as an alternative to removing gates in ApplyCircuit."
 
    Ph::usage = "Ph is the phase shift gate, which introduces phase factor exp(i*theta) upon state |1...1> of the target and control qubits. The gate is the same under different orderings of qubits, and division between control and target qubits."
 
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
 
 
 
    Begin["`Private`"]
    
    
    
        (* report a generic error that the function was passed with bad args (did not evaluate) *)
        invalidArgError[func_Symbol] := (
            Message[func::error, "Invalid arguments. See ?" <> ToString[func]];
            $Failed)
               
               
               
        (* opcodes *)
        getOpCode[gate_] :=
	        gate /. {H->0,X->1,Y->2,Z->3,Rx->4,Ry->5,Rz->6,R->7,S->8,T->9,U->10,Deph->11,Depol->12,Damp->13,SWAP->14,M->15,P->16,Kraus->17,G->18,Id->19,Ph->20,KrausNonTP->21,_->-1}
        
        (* convert MMA matrix to a flat format which can be embedded in the circuit param list *)
        codifyMatrix[matr_] :=
            Riffle[Re @ N @ Flatten @ matr, Im @ N @ Flatten @ matr]
            
        (* convert multiple MMA matrices into {#matrices, ... flattened matrices ...} *)
        codifyMatrices[matrs_] :=
            Prepend[Join @@ (codifyMatrix /@ matrs), Length @ matrs]
            
        (* sub-gate patterns *)
        pattPauli = Subscript[(X|Y|Z), _Integer];
        
        (* recognising and codifying gates into {opcode, ctrls, targs, params} *)
        gatePatterns = {
            Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[U,  (targs:__Integer)|{targs:__Integer}][matr:_List]] :> 
                {getOpCode[U], {ctrls}, {targs}, codifyMatrix[matr]},
        	Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]] :> 
                {getOpCode[gate], {ctrls}, {targs}, {args}},
        	Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]] :> 
                {getOpCode[gate], {ctrls}, {targs}, {}},
            Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][R[param_, ({paulis:pattPauli..}|Verbatim[Times][paulis:pattPauli..]|paulis:pattPauli__)]] :>
                {getOpCode[R], {ctrls}, {paulis}[[All,2]], Join[{param}, getOpCode /@ {paulis}[[All,1]]]},
            R[param_, ({paulis:pattPauli..}|Verbatim[Times][paulis:pattPauli..]|paulis:pattPauli__)] :>
                {getOpCode[R], {}, {paulis}[[All,2]], Join[{param}, getOpCode /@ {paulis}[[All,1]]]},
        	Subscript[U, (targs:__Integer)|{targs:__Integer}][matr:_List] :> 
                {getOpCode[U], {}, {targs}, codifyMatrix[matr]},
            Subscript[Kraus, (targs:__Integer)|{targs:__Integer}][matrs_List] :>
                {getOpCode[Kraus], {}, {targs}, codifyMatrices[matrs]},
            Subscript[KrausNonTP, (targs:__Integer)|{targs:__Integer}][matrs_List] :>
                {getOpCode[KrausNonTP], {}, {targs}, codifyMatrices[matrs]},
            Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__] :> 
                {getOpCode[gate], {}, {targs}, {args}},
        	Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}] :> 
                {getOpCode[gate], {}, {targs}, {}},
            G[arg_] :> 
                {getOpCode[G], {}, {}, {arg}}
        };

        (* converting gate sequence to code lists: {opcodes, ctrls, targs, params} *)
        codifyCircuit[circuit_List] :=
        	circuit /. gatePatterns // Transpose
        codifyCircuit[circuit_] :=
            codifyCircuit @ {circuit}
            
        (* checking circuit format *)
        isGateFormat[Subscript[_Symbol, (__Integer)|{__Integer}]] := True
        isGateFormat[Subscript[_Symbol, (__Integer)|{__Integer}][__]] := True
        isGateFormat[R[_, (pattPauli|{pattPauli..}|Verbatim[Times][pattPauli..])]] := True
        isGateFormat[G[_]] := True
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
                circuitProgressVar = 0;
                ApplyCircuitInternal[qureg, withBackup, showProgress, circCodes],
                ProgressIndicator[circuitProgressVar]
            ]
        ApplyCircuit[qureg_Integer, {}, OptionsPattern[ApplyCircuit]] :=
            {}
        ApplyCircuit[qureg_Integer, circuit_?isCircuitFormat, OptionsPattern[ApplyCircuit]] :=
        	With[
        		{codes = codifyCircuit[circuit]},
        		Which[
        			Not @ AllTrue[codes[[4]], NumericQ, 2],
                    Message[ApplyCircuit::error, "Circuit contains non-numerical parameters!"]; $Failed,
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
        
        (* apply the derivatives of a circuit on an initial state, storing the ersults in the given quregs *)
        extractUnitaryMatrix[Subscript[U, __Integer][u_List]] := u
        extractUnitaryMatrix[Subscript[C, (__Integer|{__Integer})][Subscript[U, __Integer][u_List]]] := u
        calcUnitaryDeriv[{param_, gate_}] := 
            D[extractUnitaryMatrix[gate], param]
        CalcQuregDerivs[circuit_?isCircuitFormat, initQureg_Integer, varVals:{(_ -> _?NumericQ) ..}, derivQuregs:{__Integer}] :=
            With[
                {varOpInds = DeleteDuplicates /@ (Position[circuit, _?(MemberQ[#])][[All, 1]]& /@ varVals[[All,1]]),
                codes = codifyCircuit[(circuit /. varVals)]}, 
                Which[
                    AnyTrue[varOpInds, Length[#]<1&],
                    Message[CalcQuregDerivs::error, "One or more variables were not present in the circuit!"]; $Failed,
                    AnyTrue[varOpInds, Length[#]>1&],
                    Message[CalcQuregDerivs::error, "One or more variables appeared multiple times in the circuit!"]; $Failed,
                    Not @ AllTrue[codes[[4]], NumericQ, 2],
                    Message[CalcQuregDerivs::error, "The circuit contained variables not assigned values in varVals!"]; $Failed,
                    True,
                    With[{unitaryGates = Select[
                        Flatten[{varVals[[All,1]], circuit[[varOpInds[[All,1]]]]}, {{2},{1}}], Not[FreeQ[#, U]] &]},
                        CalcQuregDerivsInternal[
                            initQureg, derivQuregs, Flatten[varOpInds]-1,  (* maps indices from MMA to C *)
                            unpackEncodedCircuit @ codes,
                            Flatten[codifyMatrix /@ (calcUnitaryDeriv /@ unitaryGates /. varVals)]
                        ]
                    ]
                ]
            ]
        (* error for bad args *)
        CalcQuregDerivs[___] := invalidArgError[CalcQuregDerivs]
            
        (* compute a matrix of inner products; this is used in tandem with CalcQuregDerivs to populate the Li matrix *)
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
                CalcDensityInnerProductsMatrixInternal[quregIds],
                {Length @ quregIds, Length @ quregIds}
            ]
        (* compute a real vector of density innere products *)
        CalcDensityInnerProducts[rhoId_Integer, omegaIds:{__Integer}] :=
            CalcDensityInnerProductsVectorInternal[rhoId, omegaIds]
        (* error for bad args *)
        CalcDensityInnerProducts[___] := invalidArgError[CalcDensityInnerProducts]

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
            
        (* compute the expected value of a Pauli product *)
        CalcExpecPauliProd[qureg_Integer, Verbatim[Times][paulis:pattPauli..], workspace_Integer] :=
            CalcExpecPauliProdInternal[qureg, workspace, getOpCode /@ {paulis}[[All,1]], {paulis}[[All,2]]]
        CalcExpecPauliProd[qureg_Integer, Subscript[pauli:(X|Y|Z),targ:_Integer], workspace_Integer] :=
            CalcExpecPauliProdInternal[qureg, workspace, getOpCode /@ {pauli}, {targ}]
        CalcExpecPauliProd[___] := invalidArgError[CalcExpecPauliProd]
            
        (* compute the expected value of a weighted sum of Pauli products *)
        pattXYZI = Subscript[X|Y|Z|Id, _Integer];
        getPauliSumTermCoeff[pauli:pattXYZI] = 1;
        getPauliSumTermCoeff[Verbatim[Times][coeff:_?NumericQ:1, ___]] := coeff
        getPauliSumTermCodes[pauli:pattXYZI] := {getOpCode @ pauli[[1]]}
        getPauliSumTermCodes[Verbatim[Times][___?NumericQ, paulis:pattXYZI..]] := getOpCode /@ {paulis}[[All, 1]]
        getPauliSumTermTargs[pauli:pattXYZI] := {pauli[[2]]}
        getPauliSumTermTargs[Verbatim[Times][___?NumericQ, paulis:pattXYZI ..]] := {paulis}[[All, 2]]
        (* sum of individual paulis or weighted pauli products *)        
        pattPauliSum = Verbatim[Plus][ ( pattXYZI | Verbatim[Times][___?NumericQ, pattXYZI..] ) .. ]
        
        CalcExpecPauliSum[qureg_Integer, paulis:pattPauliSum, workspace_Integer] := 
            With[{
                coeffs = getPauliSumTermCoeff /@ List @@ paulis,
                codes = getPauliSumTermCodes /@ List @@ paulis,
                targs = getPauliSumTermTargs /@ List @@ paulis
                },
                If[
                    And @@ DuplicateFreeQ /@ targs,
                    CalcExpecPauliSumInternal[qureg, workspace, coeffs, Flatten[codes], Flatten[targs], Length /@ targs],
                    (Message[CalcExpecPauliSum::error, "Pauli operators within a product must target unique qubits."]; $Failed)
                ]
            ]
        (* single term: single Pauli *)
        CalcExpecPauliSum[qureg_Integer, pauli:pattPauli, workspace_Integer] :=
            CalcExpecPauliSumInternal[qureg, workspace, {1}, {getOpCode @ pauli[[1]]}, {pauli[[2]]}, {1}]
        (* single term: pauli product, with or without coeff *)
        CalcExpecPauliSum[qureg_Integer, Verbatim[Times][coeff:_?NumericQ:1, paulis:pattPauli..], workspace_Integer] :=
            If[
                DuplicateFreeQ @ {paulis}[[All,2]],
                CalcExpecPauliSumInternal[qureg, workspace, {coeff}, getOpCode /@ {paulis}[[All,1]], {paulis}[[All,2]], {Length @ {paulis}}],
                (Message[CalcExpecPauliSum::error, "Pauli operators within a product must target unique qubits."]; $Failed)
            ]
        (* constant plus pauli sum *)
        pattConstPlusPauliSum = Verbatim[Plus][const_?NumericQ, pauliTerms:(pattPauli | Verbatim[Times][___?NumericQ, pattPauli..])..];
        CalcExpecPauliSum[qureg_Integer, blank:pattConstPlusPauliSum, workspace_Integer] := 
            (Message[CalcExpecPauliSum::error, "The Pauli sum contains a scalar. Perhaps you meant to multiply it onto an identity (Id) operator."]; $Failed)
        CalcExpecPauliSum[___] := invalidArgError[CalcExpecPauliSum]
            
        (* apply a weighted sum of Pauli products to a qureg *)
        ApplyPauliSum[inQureg_Integer, paulis:pattPauliSum, outQureg_Integer] :=
            With[{
                coeffs = getPauliSumTermCoeff /@ List @@ paulis,
                codes = getPauliSumTermCodes /@ List @@ paulis,
                targs = getPauliSumTermTargs /@ List @@ paulis
                },
                If[
                    And @@ DuplicateFreeQ /@ targs,
                    ApplyPauliSumInternal[inQureg, outQureg, coeffs, Flatten[codes], Flatten[targs], Length /@ targs],
                    (Message[ApplyPauliSum::error, "Pauli operators within a product must target unique qubits."]; $Failed)
                ]
            ]
        ApplyPauliSum[inQureg_Integer, pauli:pattPauli, outQureg_Integer] :=
            ApplyPauliSumInternal[inQureg, outQureg, {1}, {getOpCode @ pauli[[1]]}, {pauli[[2]]}, {1}]
        (* single term: pauli product, with or without coeff *)
        ApplyPauliSum[inQureg_Integer, Verbatim[Times][coeff:_?NumericQ:1, paulis:pattPauli..], outQureg_Integer] :=
            If[
                DuplicateFreeQ @ {paulis}[[All,2]],
                ApplyPauliSumInternal[inQureg, outQureg, {coeff}, getOpCode /@ {paulis}[[All,1]], {paulis}[[All,2]], {Length @ {paulis}}],
                (Message[ApplyPauliSum::error, "Pauli operators within a product must target unique qubits."]; $Failed)
            ]
        (* constant plus pauli sum *)
        ApplyPauliSum[inQureg_Integer, blank:pattConstPlusPauliSum, outQureg_Integer] := 
            (Message[ApplyPauliSum::error, "The Pauli sum contains a scalar. Perhaps you meant to multiply it onto an identity (Id) operator."]; $Failed)
        ApplyPauliSum[___] := invalidArgError[ApplyPauliSum]
        
        (* convert a symbolic expression of Pauli products into an analytic matrix *)
        getFullHilbertPauliMatrix[numQ_][Subscript[s_,q_]] := Module[
        	{m=ConstantArray[IdentityMatrix[2], numQ]},
        	m[[q+1]] = PauliMatrix[s /. {Id->0, X->1,Y->2,Z->3}];
        	If[Length[m]>1, KroneckerProduct @@ (Reverse @ m), First @ m]]
            
        SetAttributes[CalcPauliExpressionMatrix, HoldAll]
        CalcPauliExpressionMatrix[h_] := With[
        	{pauliPatt = Subscript[Id|X|Y|Z,_Integer]},
        	{hFlat = SimplifyPaulis[h]},
        	{nQb = Max[1 + Cases[{hFlat}, Subscript[(Id|X|Y|Z), q_]:>q, Infinity]]},
        	ReleaseHold[
        		HoldForm[hFlat] /. Verbatim[Times][a___, b:pauliPatt, c___] :>
        			RuleCondition @ Times[
        				Sequence @@ Cases[{a,b,c}, Except[pauliPatt]],
        				Dot @@ getFullHilbertPauliMatrix[nQb] /@ Cases[{a,b,c}, pauliPatt]
        			] /. p:pauliPatt :> RuleCondition @ getFullHilbertPauliMatrix[nQb][p]]]
        CalcPauliExpressionMatrix[___] := invalidArgError[CalcPauliExpressionMatrix]
                
        (* convert a real-weighted sum of Pauli products into a numerical matrix *)
        CalcPauliSumMatrix[paulis:pattPauliSum] := 
            With[{
                arrs=With[{
                    coeffs = getPauliSumTermCoeff /@ List @@ paulis,
                    codes = getPauliSumTermCodes /@ List @@ paulis,
                    targs = getPauliSumTermTargs /@ List @@ paulis
                    },
                    If[
                        And @@ DuplicateFreeQ /@ targs,
                        CalcPauliSumMatrixInternal[1+Max@Flatten@targs, coeffs, Flatten[codes], Flatten[targs], Length /@ targs],
                        (Message[CalcPauliSumMatrix::error, "Pauli operators within a product must target unique qubits."]; $Failed)
                    ]
                ]},
                If[
                    arrs === $Failed, arrs,
                    (#[[1]] + I #[[2]])& /@ Partition[arrs,2] // Transpose
                ]
            ]
        CalcPauliSumMatrix[blank:pattConstPlusPauliSum] := 
            (Message[CalcPauliSumMatrix::error, "The Pauli sum contains a scalar. Perhaps you meant to multiply it onto an identity (Id) operator."]; $Failed)
        CalcPauliSumMatrix[___] := invalidArgError[CalcPauliSumMatrix]
        
        (* convert a list of Pauli coefficients and codes into a weighted (symbolic) sum of products *)
        GetPauliSumFromCoeffs[addr_String] :=
            Plus @@ (#[[1]] If[ 
                    AllTrue[ #[[2;;]], PossibleZeroQ ],
                    Subscript[Id, 0],
                    Times @@ MapThread[
                    (   Subscript[Switch[#2, 0, Id, 1, X, 2, Y, 3, Z], #1 - 1] /. 
                        Subscript[Id, _] ->  Sequence[] & ), 
                        {Range @ Length @ #[[2 ;;]], #[[2 ;;]]}
                    ]
                ] &) /@ ReadList[addr, Number, RecordLists -> True];
        GetPauliSumFromCoeffs[___] := invalidArgError[GetPauliSumFromCoeffs]
        
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
                    
        CreateLocalQuESTEnv[fn_:"quest_link"] := If[
            FileExistsQ[fn], 
            Install[fn],  
            Message[CreateLocalQuESTEnv::error, "Local quest_link executable not found!"]; $Failed]
        CreateLocalQuESTEnv[__] := invalidArgError[CreateLocalQuESTEnv] (* no args is valid *)
            
        getExecFn["MacOS"|"MacOSX"] = "macos_quest_link";
        getExecFn["Windows"] = "windows_quest_link.exe";
        getExecFn["Linux"|"Unix"] = "linux_quest_link";
        CreateDownloadedQuESTEnv[os:("MacOS"|"MacOSX"|"Windows"|"Linux"|"Unix")] := 
            Module[{url,exec},
                url = "https://github.com/QTechTheory/QuESTlink/raw/master/Binaries/" <> getExecFn[os];
                exec = URLDownload[url, "quest_link"];
                If[os != "Windows", Run["chmod +x quest_link"]];
                Install[exec]
            ]
        CreateDownloadedQuESTEnv[] :=
            CreateDownloadedQuESTEnv[$OperatingSystem]
        CreateDownloadedQuESTEnv[__] := (
            Message[CreateDownloadedQuESTEnv::error, "Supported operating systems are Windows, Linux, Unix, MacOS, MacOSX."]; 
            $Failed)
                    
        DestroyQuESTEnv[link_] := Uninstall @ link
        DestroyQuESTEnv[___] := invalidArgError[DestroyQuESTEnv]
        
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
        factorPaulis[s_Plus] := Simplify /@ Plus @@@ GatherBy[List @@ s, getPauliSig] // Total
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
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][ R[arg_, Verbatim[Times][paulis:Subscript[(X|Y|Z), _Integer]..]] ]] := {Join[{R}, {paulis}[[All,1]]], {ctrls}, {paulis}[[All,2]]}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][ R[arg_, Subscript[pauli:(X|Y|Z), targ_Integer]] ]] := {{R,pauli}, {ctrls}, {targ}}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]] := {gate, {},{targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]] := {gate, {}, {targs}}
        getSymbCtrlsTargs[R[arg_, Verbatim[Times][paulis:Subscript[(X|Y|Z), _Integer]..]]] := {Join[{R}, {paulis}[[All,1]]], {}, {paulis}[[All,2]]}
        getSymbCtrlsTargs[R[arg_, Subscript[pauli:(X|Y|Z), targ_Integer]]] := {{R,pauli}, {}, {targ}}
            (* little hack to enable G[x] in GetCircuitColumns *)
            getSymbCtrlsTargs[G[x_]] := {G, {}, {}}

        (* deciding how to handle gate placement *)
        getQubitInterval[{ctrls___}, {targs___}] :=
        	Interval @ {Min[ctrls,targs],Max[ctrls,targs]}
        getNumQubitsInCircuit[circ_List] :=
        	Max[1 + Cases[{circ}, Subscript[gate_, inds__]-> Max[inds], \[Infinity]],    
        		1 + Cases[{circ}, Subscript[gate_, inds__][___] -> Max[inds], \[Infinity]]]
        isContiguousBlockGate[(SWAP|M|Rz|Ph|X|R|{R, (X|Y|Z)..})] := False
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
        drawGate[{R, rots:(X|Y|Z)..}, {}, targs_List, col_] := {
            Line[{{col+.5,Min[targs]+.5},{col+.5,Max[targs]+.5}}],
            Sequence @@ MapThread[drawGate[#1/.{X->Rx,Y->Ry,Z->Rz}, {}, {#2}, col]&, {{rots}, targs}]}
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
                    gates = compactCirc[compactFlag][subcirc /. G[_] -> Nothing];
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
        getAnalGateMatrix[Subscript[U, __][m_]] = m;
        getAnalGateMatrix[Subscript[Ph, t__][a_]] = DiagonalMatrix[ Append[ConstantArray[1, 2^Length[{t}] - 1], Exp[I a]] ];
        getAnalGateMatrix[G[a_]] := Exp[I a] {{1,0},{0,1}};
        getAnalGateMatrix[Subscript[Rx, _][a_]] = MatrixExp[-I a/2 PauliMatrix[1]]; (* KroneckerProduct doesn't have a one-arg identity overload?? Bah *)
        getAnalGateMatrix[Subscript[Ry, _][a_]] = MatrixExp[-I a/2 PauliMatrix[2]];
        getAnalGateMatrix[Subscript[Rz, _][a_]] = MatrixExp[-I a/2 PauliMatrix[3]];
        getAnalGateMatrix[Subscript[Rx, t__][a_]] := MatrixExp[-I a/2 KroneckerProduct @@ ConstantArray[PauliMatrix[1],Length[{t}]]];
        getAnalGateMatrix[Subscript[Ry, t__][a_]] := MatrixExp[-I a/2 KroneckerProduct @@ ConstantArray[PauliMatrix[2],Length[{t}]]];
        getAnalGateMatrix[Subscript[Rz, t__][a_]] := MatrixExp[-I a/2 KroneckerProduct @@ ConstantArray[PauliMatrix[3],Length[{t}]]];
        getAnalGateMatrix[R[a_, pauli_]] := MatrixExp[-I a/2 getAnalGateMatrix @ pauli];
        getAnalGateMatrix[R[a_, paulis_Times]] := MatrixExp[-I a/2 * KroneckerProduct @@ (getAnalGateMatrix /@ List @@ paulis)]
        getAnalGateMatrix[Subscript[C, __][g_]] := getAnalGateMatrix[g]
        
        (* extract ctrls from gate symbols *)
        getAnalGateControls[Subscript[C, c_List][___]] := c
        getAnalGateControls[Subscript[C, c__][___]] := {c}
        getAnalGateControls[_] := {}
            
        (* extract targets from gate symbols *)
        getAnalGateTargets[Subscript[U, t_List][_]] := t
        getAnalGateTargets[Subscript[U, t__][_]] := {t}
        getAnalGateTargets[R[_, Subscript[_, t_]]] := {t}
        getAnalGateTargets[R[_, paulis_Times]] := getAnalGateTargets /@ List @@ paulis // Flatten // Reverse
        getAnalGateTargets[Subscript[C, __][g_]] := getAnalGateTargets[g]
        getAnalGateTargets[Subscript[_, t__]] := {t}
        getAnalGateTargets[Subscript[_, t__][_]] := {t}
        
        (* convert a symbolic circuit into an analytic matrix *)
        CalcCircuitMatrix[gates_List, numQb_Integer] := With[{
        	matrices = getAnalFullMatrix[
        		getAnalGateControls@#, 
        		getAnalGateTargets@#, 
        		getAnalGateMatrix@#, numQb
            ]& /@ gates},
        	Dot @@ Reverse @ matrices
        ]
        CalcCircuitMatrix[gates_List] :=
        	CalcCircuitMatrix[gates, 1 + Max @ Cases[gates, (Subscript[_, q__]|Subscript[_, q__][__]):> Max @ q, \[Infinity]]]
        CalcCircuitMatrix[___] := invalidArgError[CalcCircuitMatrix]
        
        
        
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
        getInverseGate[Subscript[U, q__][m_?MatrixQ]] := Subscript[U, q][ConjugateTranspose[m]]
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
        tidyMatrixGate[Subscript[U, q_Integer][m_]] := Subscript[U, q][Simplify @ m]
        tidyMatrixGate[Subscript[U, q__Integer][m_]] /; OrderedQ[{q}] := Subscript[U, q][Simplify @ m]
        tidyMatrixGate[Subscript[U, q__Integer][m_]] := 
        	With[{order=Ordering[{q}]},
        		Do[
        			If[order[[i]] =!= i, Block[{tmp}, With[
        				{q1={q}[[i]], q2={q}[[order[[i]]]]}, 
        				{s=CalcCircuitMatrix[{Subscript[SWAP, i-1,order[[i]]-1]}, Length[{q}]]},
        				Return @ tidyMatrixGate @ Subscript[U, Sequence@@((({q} /. q1->tmp) /. q2->q1) /. tmp->q2)][s . m . s]]]],
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
        		g:Subscript[U, q__Integer|{q__Integer}][m_] :> tidyMatrixGate[g],
        		Subscript[C, c__][g:Subscript[U, q__Integer|{q__Integer}][m_]] :> Subscript[C, tidyInds@c][tidyMatrixGate@g],
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
        				(* multiply matrices of adjacent unitaries *)
        				{
        					{ {a___, Subscript[U, q__][m1_], b___}, {c___, Subscript[U, q__][m2_], d___} } :> Sequence[{a,Subscript[U, q][m1 . m2//Simplify],b},{c,d}],
        					{ {a___, Subscript[C, ctrl__]@Subscript[U, q__][m1_], b___}, {c___, Subscript[C, ctrl__]@Subscript[U, q__][m2_], d___} } :> Sequence[{a,Subscript[C, ctrl]@Subscript[U, q][m1 . m2//Simplify],b},{c,d}]
        				},
        				(* merge all global phases *)
        				{
        					{ {a___, G[x_], b___, G[y_], c___} } :> {G[x+y//Simplify], a, b, c},
        					{ {a___, G[x_], b___}, infix___, {c___, G[y_], d___} } :> Sequence[{G[x+y//Simplify],a,b}, infix, {c,d}]
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
        				(* remove identity matrices (qubits are sorted) *)
        				Subscript[U, q__][m_] /; m === IdentityMatrix[2^Length[{q}]] -> Nothing,
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
        			Subscript[Ph, q_][\[Pi]/2] :> Subscript[S, q],
        			Subscript[Ph, q_][\[Pi]/4] :> Subscript[T, q]
        			
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
        	Join[s1, Reverse[s1]]]
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
        		{gates = R @@@ getTrotterTerms[terms, order, reps, time]},
        		gates /. R[_, Subscript[Id, _Integer]] :> Nothing]

    End[ ]
                                       
EndPackage[]

Needs["QuEST`Option`"]

Needs["QuEST`Gate`"]

Needs["QuEST`DeviceSpec`"]

