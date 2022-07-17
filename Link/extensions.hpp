
#ifndef EXTENSIONS_H
#define EXTENSIONS_H



void extension_addAdjointToSelf(Qureg qureg);

void extension_applyImagFactor(Qureg qureg, qreal imagFac);

void extension_mixDephasingDeriv(Qureg qureg, int targetQubit, qreal probDeriv);

void extension_mixDampingDeriv(Qureg qureg, int targ, qreal prob, qreal probDeriv);



#endif // EXTENSIONS_H