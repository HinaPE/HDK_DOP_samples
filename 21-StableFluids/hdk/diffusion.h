#ifndef STABLEFLUIDS_DIFFUSION_H
#define STABLEFLUIDS_DIFFUSION_H

#include <SIM/SIM_RawField.h>

namespace HinaPE
{
void Diffuse(SIM_RawField *TARGET, const SIM_RawField *MARKER, float factor);
}

#endif //STABLEFLUIDS_DIFFUSION_H
