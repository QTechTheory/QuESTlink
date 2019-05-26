
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
    
    CreateDownloadedQuESTEnv::usage = "CreateDownloadedQuESTEnv[] downloads a MacOS-CPU-QuEST server from quest.qtechtheory.org, gives it permission to run then locally connects to it."
                    
    DestroyQuESTEnv::usage = "DestroyQuESTEnv[link] disconnects from the QuEST link, which may be the remote Igor server, clearing some QuEST function definitions (but not those provided by the QuEST package)."
    
    AddWeightedStates::usage = "AddWeightedStates[fac1, q1, fac2, q2, facOut, qOut] modifies qureg qOut to be (facOut qOut + fac1 q1 + fac2 q2)."

    SetWeightedStates::usage = "AddWeightedStates[fac1, q1, fac2, q2, qOut] modifies qureg qOut to be (fac1 q1 + fac2 q2)."
            
    DrawCircuit::usage = "
DrawCircuit[circuit] generates a circuit diagram.
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
    PackageExport[Damp]
    Damp::usage = "Damp[prob] is 1 qubit amplitude damping with the givern decay probability."
    PackageExport[SWAP]
    SWAP::usage = "SWAP is a 2 qubit gate which swaps the state of two qubits."
    PackageExport[M]
    M::usage = "M is a desctructive measurement gate which measures the indicated qubits in the Z basis."
    PackageExport[P]
    P::usage = "
P[val] is a projector onto {0,1} such that the target qubits represent val in binary (left most target takes the least significant digit in val).
P[outcomes] is a projector onto the given {0,1} outcomes. The left most qubit is set to the left most outcome"
            
    Begin["`Private`"]
               
        (* opcodes *)
        getOpCode[gate_] :=
	        gate /. {H->0,X->1,Y->2,Z->3,Rx->4,Ry->5,Rz->6,S->7,T->8,U->9,Deph->10,Depol->11,Damp->12,SWAP->13,M->14,P->15,_->-1}
        
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
        
        CreateDownloadedQuESTEnv[] := Module[{linkfile},
            SetDirectory @ NotebookDirectory[];
            linkfile = URLDownload["https://quest.qtechtheory.org/QuESTlink_MacOS_CPU", "quest_link"];
            Run["chmod +x quest_link"];
            Install[linkfile]
        ]
            
                    
        DestroyQuESTEnv[link_] := Uninstall @ link
        
        (* Im[0.] = 0, how annoying *)
        AddWeightedStates[fac1_?NumericQ, q1_Integer, fac2_?NumericQ, q2_Integer, facOut_?NumericQ, qOut_Integer] :=
            AddWeightedStatesInternal[
                Re @ N @ fac1, N @ Im @ N @ fac1, q1,
                Re @ N @ fac2, N @ Im @ N @ fac2, q2,
                Re @ N @ facOut, N @ Im @ N @ facOut, qOut
            ]
            
        SetWeightedStates[fac1_?NumericQ, q1_Integer, fac2_?NumericQ, q2_Integer, qOut_Integer] :=
            AddWeightedStatesInternal[
                Re @ N @ fac1, N @ Im @ N @ fac1, q1,
                Re @ N @ fac2, N @ Im @ N @ fac2, q2,
                N @ 0, N @ 0, qOut
            ]
        
        
        
        
        (*
         * Below is only front-end code for generating circuit diagrams from
         * from circuit the same format circuit specification
         *)
         
        (* convert symbolic gate form to {symbol, ctrls, targets} *)
        getSymbCtrlsTargs[Subscript[C, (ctrls:_Integer..)|{ctrls:_Integer..}][Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}][args__]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[C, (ctrls:_Integer..)|{ctrls:_Integer..}][Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}]]] := {gate, {ctrls}, {targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}][args__]] := {gate, {},{targs}}
        getSymbCtrlsTargs[Subscript[gate_Symbol, (targs:_Integer..)|{targs:_Integer..}]] := {gate, {}, {targs}}

        (* deciding how to handle gate placement *)
        getQubitInterval[{ctrls___}, {targs___}] :=
        	Interval @ {Min[ctrls,targs],Max[ctrls,targs]}
        getNumQubitsInCircuit[circ_List] :=
        	Max[1 + Cases[{circ}, Subscript[gate_, inds__]-> Max[inds], \[Infinity]],    
        		1 + Cases[{circ}, Subscript[gate_, inds__][___] -> Max[inds], \[Infinity]]]
        needsSpecialSwap[(SWAP|M), _List] := False
        needsSpecialSwap[label_Symbol, targs_List] :=
        	And[Length[targs] > 1, Abs[targs[[1]] - targs[[2]]] > 1]
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
        	
        (* two-qubit gate graphics *)
        drawGate[symb:(Deph|Depol), {}, {targ1_,targ2_}, col_] := {
        	EdgeForm[Dashed],
        	drawGate[If[symb===Deph,\[Phi],\[CapitalDelta]], {}, {targ1,targ2}, col]}
        drawGate[label_Symbol, {}, {targ1_,targ2_}/;Abs[targ2-targ1]===1, col_] := {
        	drawDoubleBox[Min[targ1,targ2], col],
        	Text[SymbolName@label, {col+.5,targ1+.5+.5}]}
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
        DrawCircuit[circ_List, opts:OptionsPattern[]] :=
        Module[{numCols,qubitgraphics,gategraphics},
        	{numCols,qubitgraphics,gategraphics} = generateCircuitGraphics[circ, getNumQubitsInCircuit@circ];
        	Graphics[{
                FaceForm[White], EdgeForm[Black],
                qubitgraphics, gategraphics},
                opts,
        		ImageSize -> 30 (numCols+1)
        	]
        ]
            
    End[ ]
                                       
EndPackage[]