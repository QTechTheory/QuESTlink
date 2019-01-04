
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
        ApplyCircuit[circuit, qureg] modifies qureg by applying the circuit.
        ApplyCircuit[circuit, inQureg, outQureg] leaves inQureg unchanged, but modifies outQureg to be the result of applying the circuit to inQureg."
            
    DestroyQureg::usage = "DestroyQureg[qureg] destroys the qureg associated with the given ID or symbol."
    
    GetMatrix::usage = "GetMatrix[qureg] returns the state-vector or density matrix associated with the given qureg."
            
    SetMatrix::usage = "SetMatrix[qureg, matr] modifies qureg, overwriting its statevector or density matrix with that passed."
            
    CreateRemoteQuESTEnv::usage = "CreateRemoteQuESTEnv[] connects to the remote Igor server and defines several QuEST functions, returning a link object. This should be called once. The QuEST function defintions can be cleared with DestroyQuESTEnv[link]."
                    
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
        
    Begin["`Private`"]
            
        (* disables gate commutivity, replaces Times with Dot *)
        $Pre = 
          Function[{arg}, 
           ReleaseHold[
            Hold[arg] //.  
             Times[\[Alpha]___, 
               patt : (Longest[(Subscript[_Symbol, _Integer]|Subscript[_Symbol, _Integer][___]) ..]), \[Omega]___] :> 
              Times[\[Alpha], Dot[patt], \[Omega]]], HoldAll];
               
        (* opcodes *)
        getOpCode[gate_] :=
	        gate /. {H->0,X->1,Y->2,Z->3,Rx->4,Ry->5,Rz->6,S->7,T->8}
        
        (* recognising gates *)
        gatePatterns = {
        	Subscript[C, ctrl_Integer][Subscript[gate_Symbol, targ_Integer][arg_]] :> {getOpCode[gate], ctrl, targ, arg},
        	Subscript[C, ctrl_Integer][Subscript[gate_Symbol, targ_Integer]] :> {getOpCode[gate], ctrl, targ, 0},
        	Subscript[gate_Symbol, targ_Integer][arg_] :> {getOpCode[gate], -1, targ, arg},
        	Subscript[gate_Symbol, targ_Integer] :> {getOpCode[gate], -1, targ, 0}
        };

        (* converting gate sequence to code lists: {opcodes, ctrls, targs, params} *)
        codifyCircuit[circuit_Dot] :=
        	circuit /. Dot -> List /. gatePatterns // Reverse // Transpose
                
        (* applying a sequence of symoblic gates to a qureg. ApplyCircuitInternal provided by WSTP *)
        ApplyCircuit[circuit_Dot, qureg_Integer] :=
        	With[
        		{codes = codifyCircuit[circuit]},
        		If[
        			AllTrue[ codes[[4]], NumericQ],
        			ApplyCircuitInternal[qureg, codes[[1]], codes[[2]], codes[[3]], codes[[4]]],
        			Echo["Some non-numerical arguments were passed to the backend!", "Error: "]; $Failed
        		]
        	]
        (* apply a circuit to get an output state without changing input state. CloneQureg provided by WSTP *)
        ApplyCircuit[circuit_Dot, inQureg_Integer, outQureg_Integer] :=
        	Block[{},
        		QuEST`CloneQureg[outQureg, inQureg];
        		ApplyCircuit[circuit, outQureg]
        	]

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
        	With[{flatelems = 
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
                    
        getIgorLink[port_] :=
        	LinkConnect[
        		With[{host="@igor-gpu.materials.ox.ac.uk"},
        		ToString[port] <> host <> "," <> ToString[port+1] <> host],
        		LinkProtocol->"TCPIP"]
                    
        CreateRemoteQuESTEnv[] := Install @ getIgorLink @ 55555
                    
        DestroyQuESTEnv[link_] := Uninstall @ link
            
    End[ ]
                                       
EndPackage[]