
#ifndef CIRCUITS_H
#define CIRCUITS_H

#include <string>


/*
 * Codes for Mathematica gate symbols 
 */

#define OPCODE_H 0
#define OPCODE_X 1
#define OPCODE_Y 2
#define OPCODE_Z 3
#define OPCODE_Rx 4
#define OPCODE_Ry 5
#define OPCODE_Rz 6
#define OPCODE_R 7
#define OPCODE_S 8
#define OPCODE_T 9
#define OPCODE_U 10
#define OPCODE_Deph 11
#define OPCODE_Depol 12
#define OPCODE_Damp 13
#define OPCODE_SWAP 14
#define OPCODE_M 15
#define OPCODE_P 16
#define OPCODE_Kraus 17
#define OPCODE_G 18
#define OPCODE_Id 19
#define OPCODE_Ph 20
#define OPCODE_KrausNonTP 21
#define OPCODE_Matr 22
#define OPCODE_UNonNorm 23

#define NUM_OPCODES 24

static const std::string opcodeStrings[] = {
    "H",			// OPCODE_H 0
    "X",			// OPCODE_X 1
    "Y",			// OPCODE_Y 2
    "Z",			// OPCODE_Z 3
    "Rx",			// OPCODE_Rx 4
    "Ry",			// OPCODE_Ry 5
    "Rz",			// OPCODE_Rz 6
    "R",			// OPCODE_R 7
    "S",			// OPCODE_S 8
    "T",			// OPCODE_T 9
    "U",			// OPCODE_U 10
    "Deph",			// OPCODE_Deph 11
    "Depol",		// OPCODE_Depol 12
    "Damp",			// OPCODE_Damp 13
    "SWAP",			// OPCODE_SWAP 14
    "M",			// OPCODE_M 15
    "P",			// OPCODE_P 16
    "Kraus",		// OPCODE_Kraus 17
    "G",			// OPCODE_G 18
    "Id",			// OPCODE_Id 19
    "Ph",			// OPCODE_Ph 20
    "KrausNonTP",	// OPCODE_KrausNonTP 21
    "Matr",			// OPCODE_Matr 22
    "UNonNorm"		// OPCODE_UNonNorm 23
};



/*
 * Max number of target and control qubits which can be specified 
 * for an individual gate 
 */
#define MAX_NUM_TARGS_CTRLS 100



int* local_prepareCtrlCache(int* ctrls, int numCtrls, int addTarg);



/** A single quantum gate or decoherence operator.
 * A Gate instance does not need explicit deletion.
 */
class Gate {
    private:
        
        int opcode;
        
        /** All array attributes are actually pointers to a circuit-wide array 
         * shared by all Gate instances, and hence should not be modified, and do not 
         * need individual freeing. This is to avoid superfluous translation of 
         * lists between QuESTlink's C++ and QuEST's C APIs.
         */
        int* ctrls;         int numCtrls;
        int* targs;         int numTargs;
        qreal* params;      int numParams;
        
        std::string getOpcodeStr();
    
    public:
        
        /** Initialise the gate attributes, after object creation 
         * (since the latter will happen for many Gate instances in an array).
         */
        void init(int opcode, 
            int* ctrls,     int numCtrls, 
            int* targs,     int numTargs, 
            qreal* params,  int numParams);
        
        /** Getters.
         * Warning: ctrls, targs and params are shared mutable arrays!
         */
        int getOpcode() { return opcode; }
        int* getCtrlsAddr() { return ctrls; }
        int* getTargsAddr() { return targs; }
        qreal* getParamsAddr() { return params; }
        
        /** Returns whether the gate is (at least, intended) unitary 
         * (such as Rx, H, U, UNonNorm), or not (like M, P, Matr, Damp)
         */
        bool isUnitary();
        
        /** Returns whether the gate is pure, i.e. can be performed upon 
         * statevectors (true), else whether it can only be performed upon 
         * density matrices (false)
         */
        bool isPure();
            
        /** Returns the number of outputs that this gate produces when performed 
         * in a circuit. This is the number of elements added to the outputs array 
         * in applyTo(qureg, outputs).
         */
        int getNumOutputs();
        
        /** Perform this gate upon the given qureg, modifying qureg, ignoring 
         * any outputs. This is equivalent to applyTo(qureg, NULL). 
         * @throws if the gate details are invalid
         */
        void applyTo(Qureg qureg);
        
        /** Perform this gate upon the given qureg, capturing any outputs to the 
         * front of the given outputs array (unless it is NULL). Although cast to 
         * qreal, some outputs may be integers (like the results of measurement)
         * needing later recasting.
         * @throws if the gate details are invalid
         */
        void applyTo(Qureg qureg, qreal* outputs);
        
        /** Apply the conjugate transpose of this gate upon the qureg.
         * @throws if the gate details are invalid
         * @throws if the gate has no known conjugate transpose
         */
        void applyDaggerTo(Qureg qureg);
        
        /** Apply the conjugate transpose of this gate upon the qureg.
         * @throws if the gate details are invalid
         * @throws if the gate has no known inverse
         */
        void applyInverseTo(Qureg qureg);
        
        /** Apply the derivative of this gate upon the qureg. In simple cases, 
         * the derivative is with respect to the real scalar parameter of the 
         * gate. Generally, it is with respect to a variable which itself 
         * determines the gate parameter. derivParams should contain the necessary 
         * additional scalars to determine the full derivative.
         * @throws if the gate details are invalid
         * @throws if derivParams are invalid
         */
        void applyDerivTo(Qureg qureg, qreal* derivParams, int numDerivParams);
};



/** A sequence of Gate instances.
 * A Circuit must be later deleted or fall out of scope, in which case persistent 
 * MMA arrays (pointed to by the Gate instances) are freed. 
 * Some methods below are defined in circuits.cpp, others in decoders.cpp.
 */
class Circuit {
    private:
        
        Gate* gates;
        int numGates;
        
        /** Aggregate circuit info needed by freeMMA();
         */
        int totalNumCtrls;
        int totalNumTargs;
        int totalNumParams;

        /** Destroys the MMA arrays which supply ctrls, targs and params to 
         * the gate instances. This should only be called by the destructor.
         * This method is defined in decoders.cpp.
         */
        void freeMMA();
        
    public:
        
        /** Load Gate instances from the WSTP link, populating the Circuit 
         * attributes. Calling this before the WSTP messages are sent will cause 
         * a crash. Unlike the other methods defined in circuits.cpp, this method 
         * is defined in decoders.cpp.
         */
        void loadFromMMA();
        
        /** Returns gates[ind] (does not explicitly throw for out of bounds error).
         */ 
        Gate getGate(int ind);
        
        /** Returns the number of Gate instances.
         */
        int getNumGates();
        
        /** Returns the number of gates in the circuit which would 
         * produce a non-zero number of outputs when simulated by applyTo().
         */
        int getNumGatesWithOutputs();
        
        /** Returns the total number of outputs aggregated between all gates 
         * in the circuit. This is the number of elements written to array 
         * outputs in applyTo() below.
         */
        int getTotalNumOutputs();
        
        /** Returns whether the circuit contains only unitary operations (like
         * U, UNonNorm, Rx, etc) as opposed to non-unitaries (like Matr, P, M, Damp)
         */
        bool isUnitary();
        
        /** Returns whether the circuit contains only pure gates and can hence be 
         * performed upon  statevectors (true), else whether it can only be performed 
         * upon density matrices (false)
         */
        bool isPure();
        
        /** Modify qureg by sequentially applying every gate within the circuit, 
         * with increasing index. Array outputs is modified to have its first n 
         * elements modified to the gate outputs, where n = getTotalNumOutputs(),
         * unless outputs=NULL (in which case, outputs are discarded).
         * If showProgress = true, a front-end loading bar will display the progress 
         * of the circuit simulation (via local_updateCircuitProgress()).
         */
        void applyTo(Qureg qureg, qreal* outputs, bool showProgress);
        void applyTo(Qureg qureg);
        
        /** Apply only a contiguous subset of the circuit gates to qureg, starting 
         * at index startGateInd (inclusive) and ending with endGateInd (exclusive).
         * No outputs are recorded and progress is not shown.
         */
        void applySubTo(Qureg qureg, int startGateInd, int endGateInd);
        
        /** Apply the dagger of a contiguous subset of the circuit gates to qureg.
         * The subset starts at startGateInd (inclusvie), and ends at endGateInd 
         * (exclusive). Then, each gate in the subset in REVERSE order has its 
         * dagger applied. 
         */
        void applyDaggerSubTo(Qureg qureg, int startGateInd, int endGateInd);
        
        /** Apply the inverse of a contiguous subset of the circuit gates to qureg.
         * The subset starts at startGateInd (inclusvie), and ends at endGateInd 
         * (exclusive). Then, each gate in the subset in REVERSE order has its 
         * inverse applied. 
         */
        void applyInverseSubTo(Qureg qureg, int startGateInd, int endGateInd);
        
        /** Send the given list of outputs (which must have been produced from 
         * this circuit instance via applyCircuit()) to Mathematica, formatting 
         * them into sub-lists according to the circuit structure / outputting 
         * gate types.
         */
        void sendOutputsToMMA(qreal* outputs);
        
        /** Destructor will free the persistent Mathematica arrays accesssed by 
         * the gate instances, and the gates array.
         */
        ~Circuit();
};



#endif // CIRCUITS_H