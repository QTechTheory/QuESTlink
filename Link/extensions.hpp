
#ifndef EXTENSIONS_H
#define EXTENSIONS_H



bool extension_isHermitian(Qureg qureg);

void extension_addAdjointToSelf(Qureg qureg);

void extension_applyImagFactor(Qureg qureg, qreal imagFac);

void extension_applyRealFactor(Qureg qureg, qreal realFac);

void extension_mixDephasingDeriv(Qureg qureg, int targetQubit, qreal probDeriv);

void extension_mixTwoQubitDephasingDeriv(Qureg qureg, int t1, int t2, qreal probDeriv);

void extension_mixDepolarisingDeriv(Qureg qureg, int targ, qreal probDeriv);

void extension_mixTwoQubitDepolarisingDeriv(Qureg qureg, int t1, int t2, qreal probDeriv);

void extension_mixDampingDeriv(Qureg qureg, int targ, qreal prob, qreal probDeriv);

void extension_calcExpecPauliProdsFromClassicalShadow(
    std::vector<qreal> &prodExpecVals, long numProds,
    int* sampleBases, int* sampleOutcomes, int numQb, long numSamples,
    int* pauliCodes, int* pauliTargs, int* numPaulisPerProd,
    int numBatches
); // throws



#endif // EXTENSIONS_H