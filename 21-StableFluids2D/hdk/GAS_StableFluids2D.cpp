#include "GAS_StableFluids2D.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_FieldUtils.h>
#include <GAS/GAS_ProjectNonDivergent.h>
#include <GAS/GAS_Diffuse.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);

void GAS_StableFluids2D::initializeSubclass()
{
	SIM_Data::initializeSubclass();
}
void GAS_StableFluids2D::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
}
const SIM_DopDescription *GAS_StableFluids2D::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	ACTIVATE_GAS_DENSITY
	ACTIVATE_GAS_VELOCITY
	ACTIVATE_GAS_SOURCE
	ACTIVATE_GAS_TEMPERATURE
	ACTIVATE_GAS_COLLISION
	ACTIVATE_GAS_PRESSURE
	ACTIVATE_GAS_DIVERGENCE
	PARAMETER_FLOAT(Diffusion, 0.01)
	PARAMETER_FLOAT(ImpulseFactor, 0.1)
	PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   DOP_NAME,
								   DOP_ENGLISH,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	setGasDescription(DESC);
	return &DESC;
}

void ApplyImpulsePartial(SIM_RawField *DYES, SIM_VectorField *VEL, const SIM_RawField *SOURCE, float FACTOR, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(SOURCE->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float S0 = vit.getValue();
		float d0 = DYES->getCellValue(vit.x(), vit.y(), vit.z());
		DYES->addToCell(vit.x(), vit.y(), vit.z(), FACTOR * S0);
		VEL->addToCell(vit.x(), vit.y(), vit.z(), FACTOR * UT_Vector3{0, 0.1f, 0});
	}
}
THREADED_METHOD4(, true, ApplyImpulse, SIM_RawField *, DYES, SIM_VectorField *, VEL, const SIM_RawField *, SOURCE, float, FACTOR)

void DivergencePartial(SIM_RawField *DIV, const SIM_VectorField *VEL, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(DIV->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		float divergence = 0;
		for (int axis: {0, 1})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I face0 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir0);
			UT_Vector3I face1 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir1);
			float v0 = VEL->getField(axis)->field()->getValue(face0.x(), face0.y(), face0.z());
			float v1 = VEL->getField(axis)->field()->getValue(face1.x(), face1.y(), face1.z());
			divergence += v1 - v0;
		}
		divergence /= 2;
		DIV->setCellValue(vit.x(), vit.y(), vit.z(), divergence);
	}
}
THREADED_METHOD2(, true, Divergence, SIM_RawField *, DIV, const SIM_VectorField *, VEL)

void PressureJacobiPartial(SIM_RawField *PRS1, const SIM_RawField *PRS0, const SIM_RawField *DIV, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(PRS1->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float pl = PRS0->field()->getValue(vit.x() - 1, vit.y(), vit.z());
		float pr = PRS0->field()->getValue(vit.x() + 1, vit.y(), vit.z());
		float pb = PRS0->field()->getValue(vit.x(), vit.y() - 1, vit.z());
		float pt = PRS0->field()->getValue(vit.x(), vit.y() + 1, vit.z());
		float d = DIV->field()->getValue(vit.x(), vit.y(), vit.z());
		float p = (pl + pr + pb + pt + d) / 4.f;
		PRS1->setCellValue(vit.x(), vit.y(), vit.z(), p);
	}
}
THREADED_METHOD3(, true, PressureJacobi, SIM_RawField *, PRS1, const SIM_RawField *, PRS0, const SIM_RawField *, DIV)

void SubtractGradientPartial(SIM_RawField *VEL_AXIS, const SIM_RawField *PRS, int AXIS, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(VEL_AXIS->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	float h = PRS->getVoxelSize().x();

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float v = vit.getValue();
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		constexpr int dir0 = 0, dir1 = 1;
		UT_Vector3I cell0 = SIM::FieldUtils::faceToCellMap(cell, AXIS, dir0);
		UT_Vector3I cell1 = SIM::FieldUtils::faceToCellMap(cell, AXIS, dir1);
		float p0 = PRS->field()->getValue(cell0.x(), cell0.y(), cell0.z());
		float p1 = PRS->field()->getValue(cell1.x(), cell1.y(), cell1.z());
		v -= ((p1 - p0) / h);
		vit.setValue(v);
	}
}
THREADED_METHOD3(, true, SubtractGradient, SIM_RawField *, VEL_AXIS, const SIM_RawField *, PRS, int, AXIS)
void SolvePressureJacobi(SIM_VectorField *VEL, const SIM_RawField *DIV)
{
	SIM_RawField PRS0(*DIV), PRS1(*DIV);
	PRS0.fieldNC()->constant(0);
	PRS1.fieldNC()->constant(0);
	for (int i = 0; i < 500; ++i)
	{
		if (i % 2)
			PressureJacobi(&PRS0, &PRS1, DIV);
		else
			PressureJacobi(&PRS1, &PRS0, DIV);
	}
	SubtractGradient(VEL->getXField(), &PRS0, 0);
	SubtractGradient(VEL->getYField(), &PRS0, 1);
}

void AdvectPartial(SIM_RawField *TARGET, const SIM_RawField *ORIGIN, const SIM_VectorField *FLOW, float dt, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3 pos, vel;
		ORIGIN->indexToPos(vit.x(), vit.y(), vit.z(), pos);

		vel = FLOW->getValue(pos);
		pos = pos - 0.5f * dt * vel;
		vel = FLOW->getValue(pos);
		pos = pos - 0.5f * dt * vel;

		vit.setValue(ORIGIN->getValue(pos));
	}
}
THREADED_METHOD4(, true, Advect, SIM_RawField *, TARGET, const SIM_RawField *, ORIGIN, const SIM_VectorField *, FLOW, float, dt)
void Advect(SIM_RawField *TARGET, const SIM_VectorField *FLOW, float dt)
{
	const SIM_RawField ORIGIN(*TARGET);
	Advect(TARGET, &ORIGIN, FLOW, dt);
}

bool GAS_StableFluids2D::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_SOURCE);

	if (!D || !V || !S)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	if (!D->getTwoDField() || !V->getTwoDField() || !S->getTwoDField())
	{
		addError(obj, SIM_MESSAGE, "Fields must be 2D", UT_ERROR_FATAL);
		return false;
	}

	Advect(V->getXField(), V, timestep);
	Advect(V->getYField(), V, timestep);
	Advect(D->getField(), V, timestep);
	ApplyImpulse(D->getField(), V, S->getField(), getImpulseFactor());
	SIM_RawField DIV(*D->getField());
	Divergence(&DIV, V);
	SolvePressureJacobi(V, &DIV);

	return true;
}
