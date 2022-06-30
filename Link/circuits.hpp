
#ifndef CIRCUITS_H
#define CIRCUITS_H


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




void local_applyGates(
    Qureg qureg, 
    int numOps, int* opcodes, 
    int* ctrls, int* numCtrlsPerOp, 
    int* targs, int* numTargsPerOp, 
    qreal* params, int* numParamsPerOp,
    qreal* observablesCache,
    int* finalCtrlInd, int* finalTargInd, int* finalParamInd,
    int showProgress
);



#endif // CIRCUITS_H