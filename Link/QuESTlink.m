
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
    
    ApplyCircuit::usage = "ApplyCircuit[circuit, qureg] modifies qureg by applying the circuit. Returns any measurement outcomes, grouped by M operators and ordered by their order in M.
ApplyCircuit[circuit, inQureg, outQureg] leaves inQureg unchanged, but modifies outQureg to be the result of applying the circuit to inQureg.
Accepts optional arguments WithBackup."
    ApplyCircuit::error = "`1`"
    
    CalcQuregDerivs::usage = "CalcQuregDerivs[circuit, initQureg, varVals, derivQuregs] sets the given list of (deriv)quregs to be the result of applying derivatives of the parameterised circuit to the initial state. The derivQuregs are ordered by the varVals, which should be in the format {param -> value}, where param is featured in Rx, Ry, Rz, R or U (and controlled) of the given circuit ONCE (multiple times within a U matrix is allowed). The initState is unchanged. Note Rx[theta] is allowed, but Rx[f(theta)] is not. Furthermore U matrices must contain at most one parameter."
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

    CalcPauliSumMatrix::usage = "CalcPauliSumMatrix[pauliSum] returns the matrix form of the given weighted sum of Pauli operators. The number of qubits is assumed to be the largest Pauli target."
    CalcPauliSumMatrix::error = "`1`"

    DestroyQureg::usage = "DestroyQureg[qureg] destroys the qureg associated with the given ID or symbol."
    DestroyQureg::error = "`1`"
    
    GetAmp::usage = "GetAmp[qureg, index] returns the complex amplitude of the state-vector qureg at the given index.
GetAmp[qureg, row, col] returns the complex amplitude of the density-matrix qureg at index [row, col]."
    GetAmp::error = "`1`"
    
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

    DrawCircuit::usage = "DrawCircuit[circuit] generates a circuit diagram.
DrawCircuit[circuit, numQubits] generates a circuit diagram with numQubits, useful for overriding the automated inferrence of the number of qubits if incorrect.
DrawCircuit[circuit, opts] enables Graphics options to modify the circuit diagram."
    DrawCircuit::error = "`1`"

    CalcCircuitMatrix::usage = "CalcCircuitMatrix[circuit] returns an analytic expression for the given unitary circuit, which may contain undefined symbols. The number of qubits is inferred from the circuit indices (0 to maximum specified).
CalcCircuitMatrix[circuit, numQubits] gives CalcCircuitMatrix a clue about the number of present qubits."
    CalcCircuitMatrix::error = "`1`"
    
    (*
     * optional arguments to public functions
     *)
     
    PackageExport[WithBackup]
    WithBackup::usage = "Optional argument to ApplyCircuit, indicating whether to create a backup during circuit evaluation to restore the input state in case of a circuit error. This incurs additional memory (default True). If the circuit contains no error, this option has no effect besides wasting memory."
    
    PackageExport[ShowProgress]
    ShowProgress::usage = "Optional argument to ApplyCircuit, indicating whether to show a progress bar during circuit evaluation (default False). This slows evaluation slightly."
    
    (* 
     * gate symbols, needed exporting so that their use below does not refer to a private var      
     *)

    PackageExport[H]
    H::usage = "H is the Hadamard gate."
    PackageExport[X]
    X::usage = "X is the Pauli X gate, a.k.a NOT or bit-flip gate."
    PackageExport[Y]
    Y::usage = "Y is the Pauli Y gate."
    PackageExport[Z]
    Z::usage = "Z is the Pauli Z gate."
    PackageExport[Rx]
    Rx::usage = "Rx[theta] is a rotation of theta around the x-axis of the Bloch sphere."        
    PackageExport[Ry]
    Ry::usage = "Ry[theta] is a rotation of theta around the y-axis of the Bloch sphere." 
    PackageExport[Rz]
    Rz::usage = "Rz[theta] is a rotation of theta around the z-axis of the Bloch sphere. Multiple targets enacts Exp[-i \[Theta]/2 Za ... Zc]." 
    PackageExport[R]
    R::usage = "R[theta, paulis] is the unitary Exp[-i \[Theta]/2 paulis]."   
    PackageExport[S]
    S::usage = "S is the S gate, a.k.a. PI/2 gate."
    PackageExport[T]
    T::usage = "T is the T gate, a.k.a PI/4 gate."
    PackageExport[U]
    U::usage = "U[matrix] is a general 1 or 2 qubit unitary gate, enacting the given 2x2 or 4x4 matrix."
    PackageExport[Deph]
    Deph::usage = "Deph[prob] is a 1 or 2 qubit dephasing with probability prob of error."
    PackageExport[Depol]
    Depol::usage = "Depol[prob] is a 1 or 2 qubit depolarising with probability prob of error."
    PackageExport[Damp]
    Damp::usage = "Damp[prob] is 1 qubit amplitude damping with the givern decay probability."
    PackageExport[SWAP]
    SWAP::usage = "SWAP is a 2 qubit gate which swaps the state of two qubits."
    PackageExport[M]
    M::usage = "M is a destructive measurement gate which measures the indicated qubits in the Z basis."
    PackageExport[P]
    P::usage = "P[val] is a (normalised) projector onto {0,1} such that the target qubits represent val in binary (right most target takes the least significant digit in val).
P[outcomes] is a (normalised) projector onto the given {0,1} outcomes. The left most qubit is set to the left most outcome."
    PackageExport[Kraus]
    Kraus::usage = "Kraus[ops] applies a one or two-qubit Kraus map (given as a list of Kraus operators) to a density matrix."
    PackageExport[G]
    G::usage = "G[phi] applies a global phase rotation of phi, by premultiplying Exp[i phi]."
    PcakageExport[Id]
    Id::usage = "Id is an identity gate which effects no change, but can be used for forcing gate alignment in DrawCircuit."
 
    Begin["`Private`"]
    
    
    
        (* report a generic error that the function was passed with bad args (did not evaluate) *)
        invalidArgError[func_Symbol] := (
            Message[func::error, "Invalid arguments. See ?" <> ToString[func]];
            $Failed)
               
               
               
        (* opcodes *)
        getOpCode[gate_] :=
	        gate /. {H->0,X->1,Y->2,Z->3,Rx->4,Ry->5,Rz->6,R->7,S->8,T->9,U->10,Deph->11,Depol->12,Damp->13,SWAP->14,M->15,P->16,Kraus->17,G->18,_->-1}
        
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
            R[param_, ({paulis:pattPauli..}|Verbatim[Times][paulis:pattPauli..]|paulis:pattPauli__)] :>
                {getOpCode[R], {}, {paulis}[[All,2]], Join[{param}, getOpCode /@ {paulis}[[All,1]]]},
        	Subscript[U, (targs:__Integer)|{targs:__Integer}][matr:_List] :> 
                {getOpCode[U], {}, {targs}, codifyMatrix[matr]},
            Subscript[Kraus, (targs:__Integer)|{targs:__Integer}][matrs_List] :>
                {getOpCode[Kraus], {}, {targs}, codifyMatrices[matrs]},
            Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__] :> 
                {getOpCode[gate], {}, {targs}, {args}},
        	Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}] :> 
                {getOpCode[gate], {}, {targs}, {}},
            G[arg_] :> 
                {getOpCode[G], {}, {}, {arg}}
        };

        (* converting gate sequence to code lists: {opcodes, ctrls, targs, params} *)
        codifyCircuit[circuit_List] :=
        	(circuit /. Subscript[Id,__] -> Sequence[]) /. gatePatterns // Transpose
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
        ApplyCircuit[circuit_?isCircuitFormat, qureg_Integer, OptionsPattern[ApplyCircuit]] :=
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
        ApplyCircuit[circuit_?isCircuitFormat, inQureg_Integer, outQureg_Integer, opts:OptionsPattern[ApplyCircuit]] :=
        	Block[{},
        		QuEST`CloneQureg[outQureg, inQureg];
        		ApplyCircuit[circuit, outQureg, opts]
        	]
        (* error for bad args *)
        ApplyCircuit[___] := invalidArgError[ApplyCircuit]
        
        (* apply the derivatives of a circuit on an initial state, storing the ersults in the given quregs *)
        extractUnitaryMatrix[Subscript[U, __Integer][u_List]] := u
        extractUnitaryMatrix[Subscript[C, __Integer][Subscript[U, __Integer][u_List]]] := u
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
        getPauliSumTermCoeff[pauli:pattPauli] = 1;
        getPauliSumTermCoeff[Verbatim[Times][coeff:_?NumericQ:1, ___]] := coeff
        getPauliSumTermCodes[pauli:pattPauli] := {getOpCode @ pauli[[1]]}
        getPauliSumTermCodes[Verbatim[Times][___?NumericQ, paulis:pattPauli..]] := getOpCode /@ {paulis}[[All, 1]]
        getPauliSumTermTargs[pauli:pattPauli] := {pauli[[2]]}
        getPauliSumTermTargs[Verbatim[Times][___?NumericQ, paulis:pattPauli ..]] := {paulis}[[All, 2]]
        (* sum of individual paulis or weighted pauli products *)
        pattPauliSum = Verbatim[Plus][ (pattPauli | Verbatim[Times][___?NumericQ, pattPauli..])..];
        CalcExpecPauliSum[qureg_Integer, paulis:pattPauliSum, workspace_Integer] := 
            With[{
                coeffs = getPauliSumTermCoeff /@ List @@ paulis,
                codes = getPauliSumTermCodes /@ List @@ paulis,
                targs = getPauliSumTermTargs /@ List @@ paulis
                },
                CalcExpecPauliSumInternal[qureg, workspace, coeffs, Flatten[codes], Flatten[targs], Length /@ targs]
            ]
        (* single term: single Pauli *)
        CalcExpecPauliSum[qureg_Integer, pauli:pattPauli, workspace_Integer] :=
            CalcExpecPauliSumInternal[qureg, workspace, {1}, {getOpCode @ pauli[[1]]}, {pauli[[2]]}, {1}]
        (* single term: pauli product, with or without coeff *)
        CalcExpecPauliSum[qureg_Integer, Verbatim[Times][coeff:_?NumericQ:1, paulis:pattPauli..], workspace_Integer] :=
            CalcExpecPauliSumInternal[qureg, workspace, {coeff}, getOpCode /@ {paulis}[[All,1]], {paulis}[[All,2]], {Length @ {paulis}}]
        (* constant plus pauli sum *)
        pattConstPlusPauliSum = Verbatim[Plus][const_?NumericQ, pauliTerms:(pattPauli | Verbatim[Times][___?NumericQ, pattPauli..])..];
        CalcExpecPauliSum[qureg_Integer, blank:pattConstPlusPauliSum, workspace_Integer] := 
            const + CalcExpecPauliSum[qureg, Plus @@ {pauliTerms}, workspace]
        CalcExpecPauliSum[___] := invalidArgError[CalcExpecPauliSum]
            
        (* apply a weighted sum of Pauli products to a qureg *)
        ApplyPauliSum[inQureg_Integer, paulis:pattPauliSum, outQureg_Integer] :=
            With[{
                coeffs = getPauliSumTermCoeff /@ List @@ paulis,
                codes = getPauliSumTermCodes /@ List @@ paulis,
                targs = getPauliSumTermTargs /@ List @@ paulis
                },
                ApplyPauliSumInternal[inQureg, outQureg, coeffs, Flatten[codes], Flatten[targs], Length /@ targs]
            ]
        ApplyPauliSum[inQureg_Integer, pauli:pattPauli, outQureg_Integer] :=
            ApplyPauliSumInternal[inQureg, outQureg, {1}, {getOpCode @ pauli[[1]]}, {pauli[[2]]}, {1}]
        (* single term: pauli product, with or without coeff *)
        ApplyPauliSum[inQureg_Integer, Verbatim[Times][coeff:_?NumericQ:1, paulis:pattPauli..], outQureg_Integer] :=
            ApplyPauliSumInternal[inQureg, outQureg, {coeff}, getOpCode /@ {paulis}[[All,1]], {paulis}[[All,2]], {Length @ {paulis}}]
        (* constant plus pauli sum *)
        ApplyPauliSum[inQureg_Integer, blank:pattConstPlusPauliSum, outQureg_Integer] := 
            With[{
                coeffs = Append[getPauliSumTermCoeff /@ {pauliTerms}, const], (* add const as new term... *)
                codes = Append[getPauliSumTermCodes /@ {pauliTerms}, {0}],      (* with Identity=0 Pauli code... *)
                targs = Append[getPauliSumTermTargs /@ {pauliTerms}, {0}]       (* on the first (or any) quqbit *)
                },
                ApplyPauliSumInternal[inQureg, outQureg, coeffs, Flatten[codes], Flatten[targs], Length /@ targs]
            ]
        ApplyPauliSum[___] := invalidArgError[ApplyPauliSum]
                
        (* convert a weighted sum of Pauli products into a matrix *)
        CalcPauliSumMatrix[paulis:pattPauliSum] := 
            With[{
                arrs=With[{
                    coeffs = getPauliSumTermCoeff /@ List @@ paulis,
                    codes = getPauliSumTermCodes /@ List @@ paulis,
                    targs = getPauliSumTermTargs /@ List @@ paulis
                    },
                    CalcPauliSumMatrixInternal[1+Max@Flatten@targs, coeffs, Flatten[codes], Flatten[targs], Length /@ targs]
                ]},
                (#[[1]] + I #[[2]])& /@ Partition[arrs,2] // Transpose
            ]
        CalcPauliSumMatrix[blank:pattConstPlusPauliSum] := 
            With[
                {matr=CalcPauliSumMatrix[Plus @@ {pauliTerms}]},
                matr + const IdentityMatrix @ Length @ matr 
            ]
        CalcPauliSumMatrix[___] := invalidArgError[CalcPauliSumMatrix]
        
        (* convert a list of Pauli coefficients and codes into a weighted (symbolic) sum of products *)
        GetPauliSumFromCoeffs[addr_String] :=
            Plus @@ (#[[1]] Times @@ MapThread[
                (   Subscript[Switch[#2, 0, Null, 1, X, 2, Y, 3, Z], #1 - 1] /. 
                    Subscript[Null, _] -> Sequence[] & ), 
                {Range @ Length @ #[[2 ;;]], #[[2 ;;]]}
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
        
        
        
        (*
         * Below is only front-end code for generating circuit diagrams from
         * from circuit the same format circuit specification
         *)
         
        (* convert symbolic gate form to {symbol, ctrls, targets} *)
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[C, (ctrls:__Integer)|{ctrls:__Integer}][Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}][args__]] := {gate, {},{targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:__Integer)|{targs:__Integer}]] := {gate, {}, {targs}}
        getSymbCtrlsTargs[R[arg_, Verbatim[Times][paulis:Subscript[(X|Y|Z), _Integer]..]]] := {Join[{R}, {paulis}[[All,1]]], {}, {paulis}[[All,2]]}
        getSymbCtrlsTargs[R[arg_, Subscript[pauli:(X|Y|Z), targ_Integer]]] := {{R,pauli}, {}, {targ}}

        (* deciding how to handle gate placement *)
        getQubitInterval[{ctrls___}, {targs___}] :=
        	Interval @ {Min[ctrls,targs],Max[ctrls,targs]}
        getNumQubitsInCircuit[circ_List] :=
        	Max[1 + Cases[{circ}, Subscript[gate_, inds__]-> Max[inds], \[Infinity]],    
        		1 + Cases[{circ}, Subscript[gate_, inds__][___] -> Max[inds], \[Infinity]]]
        needsSpecialSwap[(SWAP|M|Rz), _List] := False
        needsSpecialSwap[{R, (X|Y|Z)..}, _List] := False
        needsSpecialSwap[label_Symbol, targs_List] :=
        	And[Length[targs] === 2, Abs[targs[[1]] - targs[[2]]] > 1]
        getFixedThenBotTopSwappedQubits[{targ1_,targ2_}] :=
        	{Min[targ1,targ2],Min[targ1,targ2]+1,Max[targ1,targ2]}
            
        (* gate and qubit graphics primitives *)
        drawCross[targ_, col_] := {
        	Line[{{col+.5,targ+.5}-{.1,.1},{col+.5,targ+.5}+{.1,.1}}],
        	Line[{{col+.5,targ+.5}-{-.1,.1},{col+.5,targ+.5}+{-.1,.1}}]}
        drawControls[{ctrls__}, {targs__}, col_] := {
        	FaceForm[Black],
        	Table[Disk[{col+.5,ctrl+.5},.1],{ctrl,{ctrls}}],
        	With[{top=Max@{ctrls,targs},bot=Min@{ctrls,targs}},
        		Line[{{col+.5,bot+.5},{col+.5,top+.5}}]]}
        drawSingleBox[targ_, col_] :=
        	Rectangle[{col+.1,targ+.1}, {col+1-.1,targ+1-.1}]
        drawDoubleBox[targ_, col_] :=
        	Rectangle[{col+.1,targ+.1}, {col+1-.1,targ+2-.1}]
        drawQubitLines[qubits_List, col_] :=
        	Table[Line[{{col,qb+.5},{col+1,qb+.5}}], {qb,qubits}]
        drawSpecialSwapLine[targ1_, targ2_, col_] := {
        	Line[{{col,targ1+.5},{col+.1,targ1+.5}}],
        	Line[{{col+.1,targ1+.5},{col+.5-.1,targ2+.5}}],
        	Line[{{col+.5-.1,targ2+.5},{col+.5,targ2+.5}}]}
        drawSpecialSwap[targ1_,targ2_,col_] := {
        	drawSpecialSwapLine[targ1,targ2,col],
        	drawSpecialSwapLine[targ2,targ1,col]}
        	
        (* single qubit gate graphics *)
        drawGate[Id, {}, {targs___}, col_] :=
            {}
        drawGate[M, {}, {targs___}, col_] :=
        	Table[{
        		drawSingleBox[targ,col],
        		Circle[{col+.5,targ+.5-.4}, .4, {.7,\[Pi]-.7}],
        		Line[{{col+.5,targ+.5-.25}, {col+.5+.2,targ+.5+.3}}]
        		}, {targ, {targs}}]
        drawGate[Deph, {}, {targ_}, col_] := {
        	EdgeForm[Dashed], drawGate[\[Phi], {}, {targ}, col]}
        drawGate[Depol, {}, {targ_}, col_] := {
            EdgeForm[Dashed], drawGate[\[CapitalDelta], {}, {targ}, col]}
        drawGate[Damp, {}, {targ_}, col_] := {
            EdgeForm[Dashed], drawGate[\[Gamma], {}, {targ}, col]}
        drawGate[Kraus, {}, {targ_}, col_] := {
            EdgeForm[Dashed], drawGate[\[Kappa], {}, {targ}, col]}
        drawGate[X, {}, {targ_}, col_] := {
            Circle[{col+.5,targ+.5},.25],
            Line[{{col+.5,targ+.5-.25},{col+.5,targ+.5+.25}}]
        }
        drawGate[label_Symbol, {}, {targ_}, col_] := {
        	drawSingleBox[targ, col],
        	Text[SymbolName@label, {col+.5,targ+.5}]}
            
        (* special gate graphics *)
        drawGate[SWAP, {}, {targs___}, col_] := {
        	(drawCross[#,col]&) /@ {targs},
        	Line[{{col+.5,.5+Min@targs},{col+.5,.5+Max@targs}}]}
        drawGate[Z, {ctrls__}, {targ_}, col_] := {
            drawControls[{ctrls,targ},{targ},col],
            Line[{{col+.5,.5+Min@ctrls},{col+.5,.5+Max@ctrls}}]}
            
        (* multi-qubit gate graphics *)
        drawGate[Rz, {}, targs_List, col_] := {
            Line[{{col+.5,Min[targs]+.5},{col+.5,Max[targs]+.5}}],
            Sequence @@ (drawGate[Rz, {}, {#1}, col]& /@ targs)}
        drawGate[{R, rots:(X|Y|Z)..}, {}, targs_List, col_] := {
            Line[{{col+.5,Min[targs]+.5},{col+.5,Max[targs]+.5}}],
            Sequence @@ MapThread[drawGate[#1/.{X->Rx,Y->Ry,Z->Rz}, {}, {#2}, col]&, {{rots}, targs}]}
        	
        (* two-qubit gate graphics *)
        drawGate[symb:(Deph|Depol), {}, {targ1_,targ2_}, col_] := {
        	EdgeForm[Dashed],
        	drawGate[If[symb===Deph,\[Phi],\[CapitalDelta]], {}, {targ1,targ2}, col]}
        drawGate[Kraus, {}, {targ1_,targ2_}, col_] := {
        	EdgeForm[Dashed],
        	drawGate[\[Kappa], {}, {targ1,targ2}, col]}
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
        drawQubitColumn[isSpecialSwapCol_, specialSwapQubits_, numQubits_, curCol_] :=
        	If[isSpecialSwapCol,
        		(* for a special column, draw all middle lines then non-special left/right nubs *)
        		With[{nonspecial=Complement[Range[0,numQubits-1], specialSwapQubits]}, {
        			drawQubitLines[Range[0,numQubits-1],curCol+.5],
        			drawQubitLines[nonspecial,curCol],
        			drawQubitLines[nonspecial,curCol+1]}],
        		(* for a non special column, draw all qubit lines *)
        		drawQubitLines[Range[0,numQubits-1],curCol]
        	]
            
        generateCircuitGraphics[gates_List, numQubits_Integer] :=
        Module[{
        	qubitgraphics,gategraphics,
        	curCol,curSymb,curCtrls,curTargs,curInterval,curIsSpecialSwap,
        	prevIntervals,prevIsSpecialSwap,prevSpecialQubits},
        	
            (* outputs *)
        	qubitgraphics = {};
        	gategraphics = {};
        	
            (* status of whether a gate can fit in the previous column *)
        	curCol = 0;
        	prevIntervals = {};
        	prevIsSpecialSwap = False;
        	prevSpecialQubits = {};
            
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
        				(* draw qubit lines for the previous column *)
        				AppendTo[qubitgraphics, 
        					drawQubitColumn[prevIsSpecialSwap, prevSpecialQubits, numQubits, curCol]];
        				
        				(* then make a new column *)
        				curCol = curCol + If[prevIsSpecialSwap,2,1]; 
        				prevIntervals = {curInterval};
        				prevIsSpecialSwap = curIsSpecialSwap;
        				prevSpecialQubits = {};
        			),
        			(* will fit *)
        			AppendTo[prevIntervals, curInterval]
        		];
        		
        		(* record whether this gate requires special swaps *)
        		If[curIsSpecialSwap, 
        			With[{qbs=getFixedThenBotTopSwappedQubits[curTargs]},
        				AppendTo[prevSpecialQubits, qbs[[2]]];
        				AppendTo[prevSpecialQubits, qbs[[3]]]]];
        	
        		(* draw gate *)
        		AppendTo[gategraphics,
        			drawGate[
        				curSymb,curCtrls,curTargs,
        				curCol + If[prevIsSpecialSwap ~And~ Not[curIsSpecialSwap], .5, 0]]];,
        		{gate,gates}
        	];
        	
        	(* perform the final round of qubit drawing *)
        	AppendTo[qubitgraphics, 
        		drawQubitColumn[prevIsSpecialSwap, prevSpecialQubits, numQubits, curCol]];
        	
            (* return *)
        	{curCol, qubitgraphics, gategraphics}
        ]
        
        (* public function to fully render a circuit *)
        DrawCircuit[circ_List, numQubits_Integer, opts:OptionsPattern[]] :=
        Module[{numCols,qubitgraphics,gategraphics},
        	{numCols,qubitgraphics,gategraphics} = generateCircuitGraphics[DeleteCases[circ,G[_]], numQubits];
        	Graphics[{
                FaceForm[White], EdgeForm[Black],
                qubitgraphics, gategraphics},
                opts,
        		ImageSize -> 30 (numCols+1),
                Method -> {"ShrinkWrap" -> True}
        	]
        ]
        DrawCircuit[circ_List, opts:OptionsPattern[]] :=
            DrawCircuit[circ, getNumQubitsInCircuit[circ], opts]
        DrawCircuit[___] := invalidArgError[DrawCircuit]
        
        
        
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
        getAnalGateMatrix[Subscript[Y, _]] = PauliMatrix[2];
        getAnalGateMatrix[Subscript[Z, _]] = PauliMatrix[3];
        getAnalGateMatrix[Subscript[S, _]] = {{1,0},{0,I}};
        getAnalGateMatrix[Subscript[T, _]] = {{1,0},{0,Exp[I \[Pi]/4]}};
        getAnalGateMatrix[Subscript[SWAP, _,_]] = {{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};
        getAnalGateMatrix[Subscript[U, __][m_]] = m;
        getAnalGateMatrix[G[a_]] := Exp[I a] {{1,0},{0,1}};
        getAnalGateMatrix[Subscript[Rx, _][a_]] = MatrixExp[-I a/2 PauliMatrix[1]];
        getAnalGateMatrix[Subscript[Ry, _][a_]] = MatrixExp[-I a/2 PauliMatrix[2]];
        getAnalGateMatrix[Subscript[Rz, _][a_]] = MatrixExp[-I a/2 PauliMatrix[3]];
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
            
    End[ ]
                                       
EndPackage[]