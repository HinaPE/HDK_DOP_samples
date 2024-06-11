#include "source.h"

#include <SIM/SIM_FieldUtils.h>

void HinaPE::SourceSolver::EmitPartial(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *SOURCE, SIM_VectorField *FLOW, const UT_Vector3 &initV, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		fpreal value = vit.getValue();
		value += SOURCE->getValue(vit.x(), vit.y(), vit.z());
		vit.setValue(value);

		for (int axis: {0, 1, 2})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			UT_Vector3I face0 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir0);
			UT_Vector3I face1 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir1);
			float v0 = FLOW->getField(axis)->field()->getValue(face0.x(), face0.y(), face0.z());
			float v1 = FLOW->getField(axis)->field()->getValue(face1.x(), face1.y(), face1.z());
			FLOW->getField(axis)->fieldNC()->setValue(face0, v0 + initV[axis]);
			FLOW->getField(axis)->fieldNC()->setValue(face1, v1 + initV[axis]);
		}
	}
}

void HinaPE::SourceSolver::BuoyancyPartial(SIM_VectorField *FLOW, const SIM_RawField *DENSITY, const SIM_RawField *TEMPERATURE, float T_avg, float dt, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(FLOW->getYField()->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3 pos;
		FLOW->indexToPos(1, vit.x(), vit.y(), vit.z(), pos);

		float d = DENSITY->getValue(pos);
		float t = TEMPERATURE->getValue(pos);

		float value = vit.getValue();
		value += dt * (d * BuoyancySmokeDensityFactor + (t - T_avg) * BuoyancySmokeTemperatureFactor);
		vit.setValue(value);
	}
}
