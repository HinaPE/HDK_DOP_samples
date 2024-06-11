#include "diffuse.h"

#include <SIM/SIM_FieldUtils.h>

static int To1DIdx(const UT_Vector3I &idx, const UT_Vector3I &res) { return idx.x() + res.x() * (idx.y() + res.y() * idx.z()); }
static UT_Vector3I To3DIdx(int idx, const UT_Vector3I &res)
{
	UT_Vector3I ret;
	ret.z() = idx / (res.x() * res.y());
	idx -= ret.z() * res.x() * res.y();
	ret.y() = idx / res.x();
	ret.x() = idx % res.x();
	return ret;
}

float DiffusionSolver::PCG(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor)
{
	exint size = ORIGIN->numVoxels();
	UT_SparseMatrixF A(size, size);
	UT_VectorF x(0, size);
	UT_VectorF b(0, size);

	BuildLHS(&A, ORIGIN, MARKER, factor);
	BuildRHS(&b, ORIGIN);
	x = b;

	A.compile();

	UT_SparseMatrixRowF AImpl;
	AImpl.buildFrom(A);
	float error = AImpl.solveConjugateGradient(x, b, nullptr);
	WriteResult(TARGET, &x);
	return error;
}
float DiffusionSolver::GaussSeidel(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor)
{
	// TODO: Implement Gauss-Seidel
	return 0;
}
float DiffusionSolver::Jacobi(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor)
{
	// TODO: Implement Jacobi
	return 0;
}
void DiffusionSolver::BuildLHSPartial(UT_SparseMatrixF *A, const UT_VoxelArrayF *ORIGIN, const UT_VoxelArrayI *MARKER, float factor, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(MARKER);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, MARKER->getVoxelRes());
		switch ((*MARKER)(cell))
		{
			case 1: // Fluid
			{
				A->addToElement(idx, idx, 1.0f);
				for (int axis: {0, 1, 2})
				{
					constexpr int dir0 = 0, dir1 = 1;
					UT_Vector3I cell0 = SIM::FieldUtils::cellToCellMap(cell, axis, dir0);
					UT_Vector3I cell1 = SIM::FieldUtils::cellToCellMap(cell, axis, dir1);
					int idx0 = To1DIdx(cell0, MARKER->getVoxelRes());
					int idx1 = To1DIdx(cell1, MARKER->getVoxelRes());
					if (cell0[axis] >= 0)
					{
						A->addToElement(idx, idx0, -factor);
						A->addToElement(idx, idx, factor);
					}
					if (cell1[axis] < MARKER->getVoxelRes()[axis])
					{
						A->addToElement(idx, idx1, -factor);
						A->addToElement(idx, idx, factor);
					}
				}
			}
				break;
			default:
				A->addToElement(idx, idx, 1.0f);
				break;
		}
	}
}
void DiffusionSolver::BuildRHSPartial(UT_VectorF *b, const UT_VoxelArrayF *ORIGIN, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(ORIGIN);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, ORIGIN->getVoxelRes());
		(*b)(idx) = vit.getValue();
	}
}
void DiffusionSolver::WriteResultPartial(UT_VoxelArrayF *TARGET, const UT_VectorF *x, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, TARGET->getVoxelRes());
		vit.setValue((*x)(idx));
	}
}
