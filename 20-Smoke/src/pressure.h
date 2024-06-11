#ifndef SMOKE_PRESSURE_H
#define SMOKE_PRESSURE_H

#include <SIM/SIM_RawField.h>
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_ThreadedAlgorithm.h>

namespace HinaPE
{
struct PoissonSolver
{
	float PCG(UT_VoxelArrayF *PRESSURE, const UT_VoxelArrayF *DIVERGENCE, const UT_VoxelArrayI *MARKER);
	float GaussSeidel(UT_VoxelArrayF *PRESSURE, const UT_VoxelArrayF *DIVERGENCE, const UT_VoxelArrayI *MARKER);
	float Jacobi(UT_VoxelArrayF *PRESSURE, const UT_VoxelArrayF *DIVERGENCE, const UT_VoxelArrayI *MARKER);

	void BuildLHSPartial(UT_SparseMatrixF *A, const UT_VoxelArrayI *MARKER, const UT_JobInfo &info);
	THREADED_METHOD2(PoissonSolver, true, BuildLHS, UT_SparseMatrixF *, A, const UT_VoxelArrayI *, MARKER);

	void BuildRHSPartial(UT_VectorF *b, const UT_VoxelArrayF *DIVERGENCE, const UT_JobInfo &info);
	THREADED_METHOD2(PoissonSolver, true, BuildRHS, UT_VectorF *, b, const UT_VoxelArrayF *, DIVERGENCE);

	void WriteResultPartial(UT_VoxelArrayF *PRESSURE, const UT_VectorF *x, const UT_JobInfo &info);
	THREADED_METHOD2(PoissonSolver, true, WriteResult, UT_VoxelArrayF *, PRESSURE, const UT_VectorF *, x);
};
}

#endif //SMOKE_PRESSURE_H
