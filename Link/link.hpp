#ifndef LINK_H
#define LINK_H

#include "QuEST.h"
#include <vector>


/*
 * Global instance of QuESTEnv, created when MMA is linked.
 */
extern QuESTEnv env;

/*
 * Collection of instantiated Quregs
 */
extern std::vector<Qureg> quregs;
extern std::vector<bool> quregIsCreated;


#endif // LINK_H