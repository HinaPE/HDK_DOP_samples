#ifndef STABLEFLUIDS_POISSON_H
#define STABLEFLUIDS_POISSON_H

#include <SIM/SIM_VectorField.h>

namespace HinaPE
{
void ProjectToNonDivergent(SIM_VectorField *FLOW, const SIM_RawField *MARKER);
}

#endif //STABLEFLUIDS_POISSON_H
