#ifndef SMOKE_ADVECT_H
#define SMOKE_ADVECT_H

#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_RawField.h>
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_ThreadedAlgorithm.h>

struct AdvectionSolver
{
	THREADED_METHOD4(AdvectionSolver, true, SemiLagrangian, UT_VoxelArrayF *, TARGET, const SIM_RawField*, ORIGIN, const SIM_VectorField *, FLOW, float, dt);
	void SemiLagrangianPartial(UT_VoxelArrayF *TARGET, const SIM_RawField *ORIGIN, const SIM_VectorField *FLOW, float dt, const UT_JobInfo &info);
};

#endif //SMOKE_ADVECT_H
