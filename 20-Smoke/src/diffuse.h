#ifndef SMOKE_DIFFUSE_H
#define SMOKE_DIFFUSE_H

#include <SIM/SIM_RawField.h>
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_ThreadedAlgorithm.h>

struct DiffusionSolver
{
	float PCG(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor);


};

#endif //SMOKE_DIFFUSE_H
