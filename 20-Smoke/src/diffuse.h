#ifndef SMOKE_DIFFUSE_H
#define SMOKE_DIFFUSE_H

#include <SIM/SIM_RawField.h>
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_ThreadedAlgorithm.h>

namespace HinaPE
{
struct DiffusionSolver
{
	float PCG(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor);
	float GaussSeidel(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor);
	float Jacobi(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor);

	THREADED_METHOD4(DiffusionSolver, true, BuildLHS, UT_SparseMatrixF *, A, const UT_VoxelArrayF *, ORIGIN, const UT_VoxelArrayI *, MARKER, float, factor);
	void BuildLHSPartial(UT_SparseMatrixF *A, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor, const UT_JobInfo &info);

	THREADED_METHOD2(DiffusionSolver, true, BuildRHS, UT_VectorF *, b, const UT_VoxelArrayF *, ORIGIN);
	void BuildRHSPartial(UT_VectorF *b, const UT_VoxelArrayF *ORIGIN, const UT_JobInfo &info);

	THREADED_METHOD2(DiffusionSolver, true, WriteResult, UT_VoxelArrayF *, TARGET, const UT_VectorF *, x);
	void WriteResultPartial(UT_VoxelArrayF *TARGET, const UT_VectorF *x, const UT_JobInfo &info);
};
}

#endif //SMOKE_DIFFUSE_H
