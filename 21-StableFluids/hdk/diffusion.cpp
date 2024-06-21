#include "diffusion.h"

#include <SIM/SIM_FieldUtils.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_SparseMatrix.h>

static int To1DIdx(const UT_Vector3I &idx, const UT_Vector3I &res) { return idx.x() + res.x() * (idx.y() + res.y() * idx.z()); }
namespace HinaPE
{
void BuildLHSPartial(UT_SparseMatrixF *A, const SIM_RawField *ORIGIN, float factor, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(ORIGIN->field());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, ORIGIN->getVoxelRes());
		A->addToElement(idx, idx, 1.0f);
		for (int axis: {0, 1, 2})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I cell0 = SIM::FieldUtils::cellToCellMap(cell, axis, dir0);
			UT_Vector3I cell1 = SIM::FieldUtils::cellToCellMap(cell, axis, dir1);
			int idx0 = To1DIdx(cell0, ORIGIN->getVoxelRes());
			int idx1 = To1DIdx(cell1, ORIGIN->getVoxelRes());
			if (cell0[axis] >= 0)
			{
				A->addToElement(idx, idx0, -factor);
				A->addToElement(idx, idx, factor);
			}
			if (cell1[axis] < ORIGIN->getVoxelRes()[axis])
			{
				A->addToElement(idx, idx1, -factor);
				A->addToElement(idx, idx, factor);
			}
		}
	}
}
THREADED_METHOD3(, true, BuildLHS, UT_SparseMatrixF *, A, const SIM_RawField *, ORIGIN, float, factor);
void BuildRHSPartial(UT_VectorF *RHS, const SIM_RawField *ORIGIN, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(ORIGIN->field());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, ORIGIN->getVoxelRes());
		(*RHS)(idx) = vit.getValue();
	}
}
THREADED_METHOD2(, true, BuildRHS, UT_VectorF *, b, const SIM_RawField *, ORIGIN)

void WriteTargetPartial(SIM_RawField *Target, const UT_VectorF *b, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(Target->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, Target->getVoxelRes());
		vit.setValue((*b)(idx));
	}
}
THREADED_METHOD2(, true, WriteTarget, SIM_RawField *, Target, const UT_VectorF *, b)
}

void HinaPE::Diffuse(SIM_RawField *TARGET, const SIM_RawField *, float factor)
{
	exint size = TARGET->field()->numVoxels();
	UT_SparseMatrixF A(size, size);
	UT_VectorF x(0, size);
	UT_VectorF b(0, size);

	SIM_RawField ORIGIN(*TARGET);
	BuildLHSNoThread(&A, &ORIGIN, factor);
	BuildRHSNoThread(&b, &ORIGIN);
	x = b;

	A.compile();

	UT_SparseMatrixRowF AImpl;
	AImpl.buildFrom(A);
	AImpl.solveConjugateGradient(x, b, nullptr);
	WriteTarget(TARGET, &x);
}
