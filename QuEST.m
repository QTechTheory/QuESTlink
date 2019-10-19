
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
    
    ApplyCircuit::usage = "ApplyCircuit[circuit, qureg] modifies qureg by applying the circuit. Returns any measurement outcomes, grouped by M operators and ordered by their order in M.
ApplyCircuit[circuit, inQureg, outQureg] leaves inQureg unchanged, but modifies outQureg to be the result of applying the circuit to inQureg."
    
    CalcQuregDerivs::usage = "CalcQuregDerivs[circuit, initQureg, varVals, derivQuregs] sets the given list of (deriv)quregs to be the result of applying derivatives of the parameterised circuit to the initial state. The derivQuregs are ordered by the varVals, which should be in the format {param -> value}, where param is featured in Rx, Ry, Rz, R or U (and controlled) of the given circuit ONCE (multiple times within a U matrix is allowed). The initState is unchanged. Note Rx[theta] is allowed, but Rx[f(theta)] is not. Furthermore U matrices must contain at most one parameter."
    
    CalcInnerProducts::usage = "CalcInnerProducts[quregIds] returns a Hermitian matrix with i-th j-th element CalcInnerProduct[quregIds[i], quregIds[j]].
CalcInnerProducts[braId, ketIds] returns a complex vector with i-th element CalcInnerProduct[braId, ketId[i]]."
    
    Circuit::usage = "Circuit[gates] converts a product of gates into a left-to-right circuit, preserving order."
    
    Operator::usage = "Operator[gates] converts a product of gates into a right-to-left circuit."
    
    CalcExpecPauliProd::usage = "CalcExpecPauliProd[qureg, paulis, workspace] evaluates the expected value of a product of Paulis. workspace must be a qureg of equal dimensions to qureg."

    CalcExpecPauliSum::usage = "CalcExpecPauliSum[qureg, pauliSum, workspace] evaluates the expected value of a weighted sum of Pauli products, of a normalised qureg. workspace must be a qureg of equal dimensions to qureg. qureg is unchanged, and workspace is modified."

    ApplyPauliSum::usage = "ApplyPauliSum[inQureg, pauliSum, outQureg] modifies outQureg to be the result of applying the weighted sum of Paulis to inQureg."

    CalcPauliSumMatrix::usage = "CalcPauliSumMatrix[pauliSum] returns the matrix form of the given weighted sum of Pauli operators. The number of qubits is assumed to be the largest Pauli target."

    DestroyQureg::usage = "DestroyQureg[qureg] destroys the qureg associated with the given ID or symbol."
    
    GetAmp::usage = "GetAmp[qureg, index] returns the complex amplitude of the state-vector qureg at the given index.
GetAmp[qureg, row, col] returns the complex amplitude of the density-matrix qureg at index [row, col]."
    
    GetQuregMatrix::usage = "GetQuregMatrix[qureg] returns the state-vector or density matrix associated with the given qureg."
            
    SetQuregMatrix::usage = "SetQuregMatrix[qureg, matr] modifies qureg, overwriting its statevector or density matrix with that passed."
            
    CreateRemoteQuESTEnv::usage = "CreateRemoteQuESTEnv[id] connects to the remote Igor server (on port 50000+id and 50100+id) and defines several QuEST functions, returning a link object. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link]."
             
    CreateLocalQuESTEnv::usage = "CreateLocalQuESTEnv[] connects to a local Mathematica backend, running single-CPU QuEST. This requires a compatible 'quest_link' executable is in the same directory as the notebook. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link]."
    
    CreateDownloadedQuESTEnv::usage = "CreateDownloadedQuESTEnv[] downloads a MacOS-CPU-QuEST server from quest.qtechtheory.org, gives it permission to run then locally connects to it. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link]."
                    
    DestroyQuESTEnv::usage = "DestroyQuESTEnv[link] disconnects from the QuEST link, which may be the remote Igor server or a loca instance, clearing some QuEST function definitions (but not those provided by the QuEST package)."

    SetWeightedQureg::usage = "SetWeightedQureg[fac1, q1, fac2, q2, facOut, qOut] modifies qureg qOut to be (facOut qOut + fac1 q1 + fac2 q2). qOut can be one of q1 an q2.
SetWeightedQureg[fac1, q1, fac2, q2, qOut] modifies qureg qOut to be (fac1 q1 + fac2 q2). qOut can be one of q1 an q2."

    DrawCircuit::usage = "DrawCircuit[circuit] generates a circuit diagram.
DrawCircuit[circuit, numQubits] generates a circuit diagram with numQubits, useful for overriding the automated inferrence of the number of qubits if incorrect.
DrawCircuit[circuit, opts] enables Graphics options to modify the circuit diagram."
            
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
    M::usage = "M is a desctructive measurement gate which measures the indicated qubits in the Z basis."
    PackageExport[P]
    P::usage = "P[val] is a (normalised) projector onto {0,1} such that the target qubits represent val in binary (right most target takes the least significant digit in val).
P[outcomes] is a (normalised) projector onto the given {0,1} outcomes. The left most qubit is set to the left most outcome."
    PackageExport[Kraus]
    Kraus::usage = "Kraus[ops] applies a one or two-qubit Kraus map (given as a list of Kraus operators) to a density matrix."
    PackageExport[G]
    G::usage = "G[phi] applies a global phase rotation of phi, by premultiplying Exp[i phi]."
 
    Begin["`Private`"]
               
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
                
        (* applying a sequence of symoblic gates to a qureg. ApplyCircuitInternal provided by WSTP *)
        ApplyCircuit[circuit_?isCircuitFormat, qureg_Integer] :=
        	With[
        		{codes = codifyCircuit[circuit]},
        		If[
        			AllTrue[codes[[4]], NumericQ, 2],
        			ApplyCircuitInternal[qureg, unpackEncodedCircuit[codes]], 
        			Echo["Circuit contains non-numerical parameters!", "Error: "]; $Failed
        		]
        	]
        (* apply a circuit to get an output state without changing input state. CloneQureg provided by WSTP *)
        ApplyCircuit[circuit_?isCircuitFormat, inQureg_Integer, outQureg_Integer] :=
        	Block[{},
        		QuEST`CloneQureg[outQureg, inQureg];
        		ApplyCircuit[circuit, outQureg]
        	]
        
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
                    Echo["One or more variables were not present in the circuit!", "Error: "]; $Failed,
                    AnyTrue[varOpInds, Length[#]>1&],
                    Echo["One or more variables appeared multiple times in the circuit!", "Error: "]; $Failed,
                    Not @ AllTrue[codes[[4]], NumericQ, 2],
                    Echo["The circuit contained variables not assigned values in varVals!", "Error: "]; $Failed,
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

        (* destroying a remote qureg, and clearing the local symbol *)
        SetAttributes[DestroyQureg, HoldAll];
        DestroyQureg[qureg_Integer] :=
        	DestroyQuregInternal[qureg]
        DestroyQureg[qureg_Symbol] :=
        	Block[{}, DestroyQuregInternal[ReleaseHold@qureg]; Clear[qureg]]

        (* get a local matrix representation of the remote qureg. GetStateVecInternal provided by WSTP *)
        GetQuregMatrix[qureg_Integer] :=
        	With[{data = GetStateVecInternal[qureg]},
        		Which[
        			data === -1,
        			$Failed,
        			data[[2]] === 0,
        			MapThread[#1 + I #2 &, {data[[3]], data[[4]]}],
        			data[[2]] === 1,
        			Transpose @ ArrayReshape[
        				MapThread[#1 + I #2 &, {data[[3]], data[[4]]}], 
        				{2^data[[1]],2^data[[1]]}]
        		]
        	]

        (* overwrite the state of a remote qureg. InitStateFromAmps provided by WSTP *)
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
            
        (* compute the expected value of a Pauli product *)
        CalcExpecPauliProd[qureg_Integer, Verbatim[Times][paulis:pattPauli..], workspace_Integer] :=
            CalcExpecPauliProdInternal[qureg, getOpCode /@ {paulis}[[All,1]], {paulis}[[All,2]], workspace]
        CalcExpecPauliProd[qureg_Integer, Subscript[pauli:(X|Y|Z),targ:_Integer], workspace_Integer] :=
            CalcExpecPauliProdInternal[qureg, getOpcode /@ {pauli}, {targ}, workspace]
            
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
        
        getIgorLink[id_] :=
        	LinkConnect[
        		With[{host="@129.67.85.74",startport=50000},
        		ToString[startport+id] <> host <> "," <> ToString[startport+id+100] <> host],
        		LinkProtocol->"TCPIP"]
                    
        CreateRemoteQuESTEnv[port_Integer] := Install @ getIgorLink @ port
                    
        CreateLocalQuESTEnv[] := With[{fn="quest_link"},
            If[FileExistsQ[fn], Install[fn],  Echo["Local quest_link executable not found!", "Error: "]; $Failed]
        ]
        
        CreateDownloadedQuESTEnv[] := Module[{linkfile},
            SetDirectory @ NotebookDirectory[];
            linkfile = URLDownload["https://quest.qtechtheory.org/QuESTlink_MacOS_CPU", "quest_link"];
            Run["chmod +x quest_link"];
            Install[linkfile]
        ]
            
                    
        DestroyQuESTEnv[link_] := Uninstall @ link
        
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
                N @ 0, N @ 0, qOut
            ]
        
        GetAmp[qureg_Integer, index_Integer] := GetAmpInternal[qureg, index, -1]
        GetAmp[qureg_Integer, row_Integer, col_Integer] := GetAmpInternal[qureg, row, col]
        
        
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
            
    End[ ]
                                       
EndPackage[]