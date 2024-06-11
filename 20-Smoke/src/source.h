#ifndef SMOKE_SOURCE_H
#define SMOKE_SOURCE_H

#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>

namespace HinaPE
{
struct SourceSolver
{
	THREADED_METHOD4(SourceSolver, true, Emit, UT_VoxelArrayF *, TARGET, const UT_VoxelArrayF *, SOURCE, SIM_VectorField *, FLOW, const UT_Vector3&, initV);
	void EmitPartial(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *SOURCE, SIM_VectorField *FLOW, const UT_Vector3 &initV, const UT_JobInfo &info);

	THREADED_METHOD5(SourceSolver, true, Buoyancy, SIM_VectorField *, FLOW, const SIM_RawField *, DENSITY, const SIM_RawField *, TEMPERATURE, float, T_avg, float, dt);
	void BuoyancyPartial(SIM_VectorField *FLOW, const SIM_RawField *DENSITY, const SIM_RawField *TEMPERATURE, float T_avg, float dt, const UT_JobInfo &info);

	float BuoyancySmokeDensityFactor = -0.000625;
	float BuoyancySmokeTemperatureFactor = 2.0;
};
}

#endif //SMOKE_SOURCE_H
