#include "poisson.h"

#include <SIM/SIM_FieldUtils.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_SparseMatrix.h>

static int To1DIdx(const UT_Vector3I &idx, const UT_Vector3I &res) { return idx.x() + res.x() * (idx.y() + res.y() * idx.z()); }
namespace HinaPE
{
void BuildLHSPartial(UT_SparseMatrixF *A, const SIM_RawField *PRS, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(PRS->field());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, PRS->getVoxelRes());
//		A->addToElement(idx, idx, 1.0f);
		for (int axis: {0, 1, 2})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I cell0 = SIM::FieldUtils::cellToCellMap(cell, axis, dir0);
			UT_Vector3I cell1 = SIM::FieldUtils::cellToCellMap(cell, axis, dir1);
			int idx0 = To1DIdx(cell0, PRS->getVoxelRes());
			int idx1 = To1DIdx(cell1, PRS->getVoxelRes());
			if (cell0[axis] >= 0)
			{
				A->addToElement(idx, idx0, -1.0f);
				A->addToElement(idx, idx, 1.0f);
			}
			if (cell1[axis] < PRS->getVoxelRes()[axis])
			{
				A->addToElement(idx, idx1, -1.0f);
				A->addToElement(idx, idx, 1.0f);
			}
		}
	}
}
THREADED_METHOD2(, true, BuildLHS, UT_SparseMatrixF *, A, const SIM_RawField *, PRS);
void BuildRHSPartial(UT_VectorF *RHS, SIM_RawField *DIV, const SIM_VectorField *FLOW, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(DIV->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	float h = DIV->getVoxelSize().x();

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, DIV->getVoxelRes());
		float divergence = 0;
		for (int axis: {0, 1, 2})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I face0 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir0);
			UT_Vector3I face1 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir1);
			float v0 = FLOW->getField(axis)->field()->getValue(face0.x(), face0.y(), face0.z());
			float v1 = FLOW->getField(axis)->field()->getValue(face1.x(), face1.y(), face1.z());
			divergence += (v1 - v0);
		}
		divergence *= h;
		vit.setValue(divergence);
		(*RHS)(idx) = -divergence;
	}
}
THREADED_METHOD3(, true, BuildRHS, UT_VectorF *, b, SIM_RawField *, DIV, const SIM_VectorField *, FLOW)

void WritePressurePartial(SIM_RawField *PRS, const UT_VectorF *b, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(PRS->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		int idx = To1DIdx(cell, PRS->getVoxelRes());
		vit.setValue((*b)(idx));
	}
}
THREADED_METHOD2(, true, WritePressure, SIM_RawField *, PRS, const UT_VectorF *, b)

void SubtractGradientPartial(SIM_RawField *FLOW, const SIM_RawField *PRS, int AXIS, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(FLOW->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	float h = PRS->getVoxelSize().x();

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I face(vit.x(), vit.y(), vit.z());
		constexpr int dir0 = 0, dir1 = 1;
		UT_Vector3I cell0 = SIM::FieldUtils::faceToCellMap(face, AXIS, dir0);
		UT_Vector3I cell1 = SIM::FieldUtils::faceToCellMap(face, AXIS, dir1);
		float p0 = PRS->field()->getValue(cell0.x(), cell0.y(), cell0.z());
		float p1 = PRS->field()->getValue(cell1.x(), cell1.y(), cell1.z());

		float v = vit.getValue();
		v -= ((p1 - p0) / h);
		vit.setValue(v);
	}
}
THREADED_METHOD3(, true, SubtractGradient, SIM_RawField *, FLOW, const SIM_RawField *, PRS, int, AXIS)
}

void HinaPE::ProjectToNonDivergent(SIM_VectorField *FLOW, const SIM_RawField *MARKER)
{
	exint size = MARKER->field()->numVoxels();
	UT_SparseMatrixF A(size, size);
	UT_VectorF x(0, size);
	UT_VectorF b(0, size);

	SIM_RawField DIV(*MARKER);
	SIM_RawField PRS(*MARKER);
	BuildLHSNoThread(&A, &PRS);
	BuildRHSNoThread(&b, &DIV, FLOW);
	x = b;

	A.compile();

	UT_SparseMatrixRowF AImpl;
	AImpl.buildFrom(A);
	AImpl.solveConjugateGradient(x, b, nullptr);
	WritePressure(&PRS, &x);
	for (int axis: {0, 1, 2})
		SubtractGradientNoThread(FLOW->getField(axis), &PRS, axis);
}
