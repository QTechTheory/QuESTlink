
#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "QuEST_complex.h"

#include "utilities.hpp"


 
/** A single term among the partial derivatives of a parameterised circuit, 
 * after expansion via the chain rule.
 */
class DerivTerm {
    private:
        
        Gate gate;
        
        /** The index of this term's corresponding gate relative to the 
         * circuit in which it is embedded.
         */
        int gateInd;
        
        /** The index of this term's corresponding differential variable. This 
         * maps one-to-one with a Qureg to be populated with a derivative state.
         */
        int varInd;
        
        /** Additional parameter data needed to effect the derivative of this 
         * term's corresponding gate. Note all DerivTerm instances secretly share
         * the same derivParams array with different subpointers, and ergo none 
         * should attempt to modify it nor delete it.
         */
        qreal* derivParams;
        
        /** The length of derivParams, used only for on-the-fly internal error
         * checking during DerivTerm::applyTo()
         */
        int numDerivParams;
        
    public:
        
        /** Initialise the DerivTerm attributes, after object creation 
         * (since the latter will happen for many DerivTerm instances in an array).
         */
        void init(Gate gate, int gateInd, int varInd, qreal* derivParams, int numDerivParams);
        
        /** Getters.
         * Warning, derivParams is a mutable array shared between DerivTerm instances!
         */
        int getGateInd() { return gateInd; };
        int getVarInd() { return varInd; };
        qreal* getDerivParamsAddr() { return derivParams; };
        
        /** Apply this derivative term (as the analytic derivative of its corresponding 
         * gate) to qureg.
         * @throws QuESTException without qureg modification if gate or derivParams 
         *      contains an invalid gate or deriv spec
         */
        void applyTo(Qureg qureg);
};



/** A complete specification of the derivatives of a parameterised circuit with
 * respect to a collection of real variables. This can include multi-variable
 * parameter gates, gates whose parameters are functions of the differential 
 * parameters (that are already numerically evaluated), circuits with variables 
 * repeated over multiple gates, gates with general element-wise parameters 
 * (like U), and even differential-parameterised decoherence channels.
 */
class DerivCircuit {
    private:
        
        /** A complete specification of the original non-differentiated circuit 
         */
        Circuit *circuit;
        
        /** The number of represented differential variables.
         */
        int numVars;
        
        /** A complete description of the differential terms of the circuit,
         * after expansion via the chain rule. The gateInd of each successive 
         * term is increasing or the same, to permit optimisations. All instances 
         * secretly share the same mutable derivParams array loaded from MMA.
         */
        DerivTerm* terms;
        int numTerms;
        
        /** The total number of derivParams aggregated between all terms of the 
         * differential circuit. This is needed to later free the MMA-loaded 
         * derivParams array which is shared between DerivTerm instances.
         */
        int totalNumDerivParams;
        
        /** Qureg type-specific implementations of public methods 
         */
        void calcDerivEnergiesStateVec(qreal* energies, PauliHamil hamil, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        void calcDerivEnergiesDensMatr(qreal* energyGrad, PauliHamil hamil, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        qmatrix calcMetricTensorStateVec(Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        qmatrix calcMetricTensorDensMatr(Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        
        /** Destroys the MMA arrays shared between DerivTerm instances (derivPArams), 
         * invoked during the destructor. This method is defined in decoders.cpp.
         */
        void freeMMA();
        
    public:
        
        /** Load attributes from the WSTP link, including the non-differentiated 
         * circuit, and each DerivTerm info. Calling this before the WSTP messages 
         * are sent will cause an unpreventable crash.
         * Unlike the other methods defined in derivatives.cpp, this method 
         * is defined in decoders.cpp.
         */
        void loadFromMMA();
        
        /** Getters 
         */
        int getNumVars() { return numVars; };
        
        /** Modify quregs to be the derivative of the circuit state produced 
         * by attribute circuit upon constant initial state initQureg. 
         * @param quregs the list of to-be-modified quregs which match the order 
         *               of varInd between the DerivTerm instances
         * @precondition varInds between all terms lie in [0, numQuregs)
         * @precondition numQuregs = number of unique varInd between terms
         * @precondition gateInd between terms is increasing (or repeating)
         * @precondition all quregs are created, and of equal dimension & type
         * @throws QuESTException if the circuit of deriv info is invalid 
         *         (i.e. contains invalid gate details), or if the user aborts
         */
        void applyTo(Qureg* quregs, int numQuregs, Qureg initQureg, Qureg workspace);
        
        /** Modifies eneryGrad to be the gradient of the expected energy under
         * the given Hamiltonian, as prescribed by the circuit derivatives.
         * This function uses the bespoke O(#parameters) time and O(1) memory 
         * algorithm from  arXiv 2009.02823, using a novel adaptation
         * for density matrices. Note that the density-matrix version involves 
         * populating one of the workQuregs with a dense representation of the 
         * hamil; use calcDerivEnergiesDenseHamil() to pre-prepare this qureg.
         * @param energyGrad must be a pre-allocated length-numVars array.
         */ 
        void calcDerivEnergies(qreal* energyGrad, PauliHamil hamil, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        
        /** This is a density-matrix only version of calcDerivEnergies(), where 
         * hamilQureg has been pre-prepared to be a matrix form of a PauliHamil,
         * via setQuregToPauliString().
         */
        void calcDerivEnergiesDenseHamil(qreal* energyGrad, Qureg hamilQureg, Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        
        /** Returns a derivative metric tensor which can be used for quantum 
         * natural gradient. If initQureg is a state-vector and the circuit is pure, 
         * then this function computes the quantum geometric tensor; a complex matrix,
         * related to the Fubini-Study metric, the classical Fisher information matrix, 
         * and the imaginary-time Li tensor plus Berry connections. If initQureg 
         * is a density matrix and/or the circuit contains decoherence channels, this 
         * function computes the Hilbert-Schmidt distance between the derivatives 
         * of the output states, HS_ij = Tr(d rho/di d rho/dj), which is a real matrix 
         * which approximates (a scaling of) the (also real) quantum Fisher information 
         * matrix in typical NISQ settings (such as when modellable as global depolarisation,
         * see https://arxiv.org/abs/1912.08660). In both scenarios, the metric tensor 
         * is computed in O(#parameters^2) time and O(1) memory, using my algorithm 
         * from https://arxiv.org/abs/2011.02991 and a density-matrix adaptation.
         * @throws exception when numWorkQuregs != 4
         */
        qmatrix calcMetricTensor(Qureg initQureg, Qureg* workQuregs, int numWorkQuregs);
        
        /** Returns the number of working registers needed to perform the method 
         * indicated by funcName upon given the initial register.
         */
        int getNumNeededWorkQuregsFor(std::string funcName, Qureg initQureg);
        
        /** Throws an exception if the workQuregIds are invalid or if there are too 
         * few for the given method.
         * @precondition initQuregId must be valid
         */
        void validateWorkQuregsFor(std::string methodName, int initQuregId, int* workQuregIds, int numWorkQuregs);
        
        /** Destructor will free the persistent Mathematica arrays accesssed by 
         * the DerivTerm instances, as well as the Circuit. 
         */
        ~DerivCircuit();
};



#endif // DERIVATIVES_H