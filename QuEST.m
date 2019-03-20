
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
    
    ApplyCircuit::usage = "
ApplyCircuit[circuit, qureg] modifies qureg by applying the circuit. Returns any measurement outcomes.
ApplyCircuit[circuit, inQureg, outQureg] leaves inQureg unchanged, but modifies outQureg to be the result of applying the circuit to inQureg."
    
    Circuit::usage = "Circuit[gates] converts a product of gates into a left-to-right circuit, preserving order."
    
    Operator::usage = "Operator[gates] converts a product of gates into a right-to-left circuit."
    
    DestroyQureg::usage = "DestroyQureg[qureg] destroys the qureg associated with the given ID or symbol."
    
    GetMatrix::usage = "GetMatrix[qureg] returns the state-vector or density matrix associated with the given qureg."
            
    SetMatrix::usage = "SetMatrix[qureg, matr] modifies qureg, overwriting its statevector or density matrix with that passed."
            
    CreateRemoteQuESTEnv::usage = "CreateRemoteQuESTEnv[id] connects to the remote Igor server (on port 50000+id and 50100+id) and defines several QuEST functions, returning a link object. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link]."
             
    CreateLocalQuESTEnv::usage = "CreateLocalQuESTEnv[] connects to a local Mathematica backend, running single-CPU QuEST."
                    
    DestroyQuESTEnv::usage = "DestroyQuESTEnv[link] disconnects from the QuEST link, which may be the remote Igor server, clearing some QuEST function definitions (but not those provided by the QuEST package)."
            
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
    Rz::usage = "Rz[theta] is a rotation of theta around the z-axis of the Bloch sphere." 
    PackageExport[S]
    S::usage = "S is the S gate, a.k.a. PI/2 gate."
    PackageExport[T]
    T::usage = "T is the T gate, a.k.a PI/4 gate."
    PackageExport[U]
    U::usage = "U[matrix] is a general single-qubit unitary gate, enacting the given 2x2 matrix."
    PackageExport[Deph]
    Deph::usage = "Deph[prob] is a 1 or 2 qubit dephasing with probability prob of error."
    PackageExport[Depol]
    Depol::usage = "Depol[prob] is a 1 or 2 qubit depolarising with probability prob of error."
    PackageExport[SWAP]
    SWAP::usage = "SWAP is a 2 qubit gate which swaps the state of two qubits."
    PackageExport[M]
    M::usage = "M is a desctructive measurement gate which measures the indicated qubits in the Z basis"
            
    Begin["`Private`"]
               
        (* opcodes *)
        getOpCode[gate_] :=
	        gate /. {H->0,X->1,Y->2,Z->3,Rx->4,Ry->5,Rz->6,S->7,T->8,U->9,Deph->10,Depol->11,SWAP->12,M->13,_->-1}
        
        (* convert MMA matrix to QuESTs ComplexMatrix2 *)
        codifyMatrix[List[List[r0c0_, r0c1_], List[r1c0_, r1c1_]]] :=
            {Re[r0c0],Im[r0c0],
             Re[r0c1],Im[r0c1],
             Re[r1c0],Im[r1c0],
             Re[r1c1],Im[r1c1]}
        
        (* recognising gates *)
        gatePatterns = {
            Subscript[C, (ctrls:_Integer..)|{ctrls:_Integer..}][Subscript[U,  (targs:_Integer..)|{targs:_Integer..}][matr:_List]] :> 
                {getOpCode[U], {ctrls}, {targs}, codifyMatrix[matr]},
        	Subscript[C, (ctrls:_Integer..)|{ctrls:_Integer..}][Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}][args__]] :> 
                {getOpCode[gate], {ctrls}, {targs}, {args}},
        	Subscript[C, (ctrls:_Integer..)|{ctrls:_Integer..}][Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}]] :> 
                {getOpCode[gate], {ctrls}, {targs}, {}},
        	Subscript[U, (targs:_Integer..)|{targs:_Integer..}][matr:_List] :> 
                {getOpCode[U], {}, {targs}, codifyMatrix[matr]},
            Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}][args__] :> 
                {getOpCode[gate], {}, {targs}, {args}},
        	Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}] :> 
                {getOpCode[gate], {}, {targs}, {}}
        };

        (* converting gate sequence to code lists: {opcodes, ctrls, targs, params} *)
        codifyCircuit[circuit_List] :=
        	circuit /. gatePatterns // Transpose
        codifyCircuit[circuit_] :=
            codifyCircuit @ {circuit}
            
        (* checking circuit format *)
        isGateFormat[Subscript[_Symbol, (_Integer..)|{_Integer..}]] := True
        isGateFormat[Subscript[_Symbol, (_Integer..)|{_Integer..}][__]] := True
        isGateFormat[___] := False
        isCircuitFormat[circ_List] := AllTrue[circ,isGateFormat]
        isCircuitFormat[circ_?isGateFormat] := True
        isCircuitFormat[___] := False
                
        (* applying a sequence of symoblic gates to a qureg. ApplyCircuitInternal provided by WSTP *)
        ApplyCircuit[circuit_?isCircuitFormat, qureg_Integer] :=
        	With[
        		{codes = codifyCircuit[circuit]},
        		If[
        			AllTrue[codes[[4]], NumericQ, 2],
        			ApplyCircuitInternal[
                        qureg, codes[[1]], 
                        Flatten @ codes[[2]], Length /@ codes[[2]], 
                        Flatten @ codes[[3]], Length /@ codes[[3]],
                        Flatten[N /@ codes[[4]]], Length /@ codes[[4]]
                    ], 
        			Echo["Circuit contains non-numerical parameters!", "Error: "]; $Failed
        		]
        	]
        (* apply a circuit to get an output state without changing input state. CloneQureg provided by WSTP *)
        ApplyCircuit[circuit_?isCircuitFormat, inQureg_Integer, outQureg_Integer] :=
        	Block[{},
        		QuEST`CloneQureg[outQureg, inQureg];
        		ApplyCircuit[circuit, outQureg]
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
        GetMatrix[qureg_Integer] :=
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
        SetMatrix[qureg_Integer, elems_List] :=
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
        
        getIgorLink[id_] :=
        	LinkConnect[
        		With[{host="@129.67.85.74",startport=50000},
        		ToString[startport+id] <> host <> "," <> ToString[startport+id+100] <> host],
        		LinkProtocol->"TCPIP"]
                    
        CreateRemoteQuESTEnv[port_Integer] := Install @ getIgorLink @ port
                    
        CreateLocalQuESTEnv[] := Install[NotebookDirectory[] <> "quest_link"]
                    
        DestroyQuESTEnv[link_] := Uninstall @ link
            
    End[ ]
                                       
EndPackage[]