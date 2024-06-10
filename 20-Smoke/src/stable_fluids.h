#ifndef SMOKE_STABLE_FLUIDS_H
#define SMOKE_STABLE_FLUIDS_H

#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>

struct StableFluids
{
	void diffusionPartial(SIM_RawField *Density, const SIM_VectorField *Velocity, const UT_JobInfo &info);
	THREADED_METHOD2(StableFluids, Density->shouldMultiThread(), diffusion, SIM_RawField *, Density, const SIM_VectorField *, Velocity);
};

#endif //SMOKE_STABLE_FLUIDS_H
