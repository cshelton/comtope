#ifndef ENSURECTBN_H
#define ENSURECTBN_H

#include "ctbn.h"
#include "bn.h"
#include "markov.h"
#include "multirv.h"
#include "ctbn.h"

namespace ctbn {
ENSURECLASS(CTBN)
ENSURECLASS(BN)
ENSURECLASS(Markov)
ENSURECLASS(CTBNDyn)
ENSURECLASS(MultiRV)
ENSURECLASS(MarkovDyn)
ENSURECLASS(MarkovSimple)
ENSURECLASS(MarkovSimpleToggle)
ENSURECLASS1(RVCondSimpleComp,MultiZSimple)
}

#endif
