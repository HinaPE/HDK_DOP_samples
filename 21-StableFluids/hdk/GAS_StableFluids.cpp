#include "GAS_StableFluids.h"

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
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_SparseMatrix.h>

#include "poisson.h"
#include "diffusion.h"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);

#define GAS_NAME_VELOCITY_SWAP        "velocity_swap"
#define ACTIVATE_GAS_VELOCITY_SWAP static PRM_Name VelocitySwapName(GAS_NAME_VELOCITY_SWAP, "VelocitySwap"); static PRM_Default VelocitySwapNameDefault(0, GAS_NAME_VELOCITY_SWAP); PRMs.emplace_back(PRM_STRING, 1, &VelocitySwapName, &VelocitySwapNameDefault);

void GAS_StableFluids::initializeSubclass() { SIM_Data::initializeSubclass(); }
void GAS_StableFluids::makeEqualSubclass(const SIM_Data *source) { SIM_Data::makeEqualSubclass(source); }
const SIM_DopDescription *GAS_StableFluids::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	ACTIVATE_GAS_DENSITY
	ACTIVATE_GAS_VELOCITY
	ACTIVATE_GAS_VELOCITY_SWAP
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

void EmitSourcePartial(SIM_RawField *TARGET, const SIM_RawField *SOURCE, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float value = SOURCE->field()->getValue(vit.x(), vit.y(), vit.z());
		float old = vit.getValue();
		if (value > 0.1f)
			vit.setValue(value);
	}
}
THREADED_METHOD2(, true, EmitSource, SIM_RawField *, TARGET, const SIM_RawField *, SOURCE)

constexpr float BuoyancySmokeDensityFactor = -0.000625;
constexpr float BuoyancySmokeTemperatureFactor = 2.0;
void BuoyancyPartial(SIM_RawField *V, const SIM_RawField *DENSITY, const SIM_RawField *TEMPERATURE, float T_avg, float dt, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(V->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3 pos;
		V->indexToPos(vit.x(), vit.y(), vit.z(), pos);

		float d = DENSITY->getValue(pos);
		float t = TEMPERATURE->getValue(pos);

		float value = vit.getValue();
		value += dt * (d * BuoyancySmokeDensityFactor + (t - T_avg) * BuoyancySmokeTemperatureFactor);
		vit.setValue(value);
	}
}
THREADED_METHOD5(, true, Buoyancy, SIM_RawField *, FLOW, const SIM_RawField *, DENSITY, const SIM_RawField *, TEMPERATURE, float, T_avg, float, dt);

bool GAS_StableFluids::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_ScalarField *T = getScalarField(obj, GAS_NAME_TEMPERATURE);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_SOURCE);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);

	if (!D || !T || !V || !S)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}
	SIM_RawField MARKER(*D->getField());
	MARKER.fieldNC()->constant(0);
	float h = D->getField()->getVoxelSize().x();
	float factor = getDiffusion() * timestep / (h * h);


	// Velocity Step
	Buoyancy(V->getYField(), D->getField(), T->getField(), T->getField()->average(), timestep);
//	for (int axis: {0, 1, 2})
//		HinaPE::Diffuse(V->getField(axis), &MARKER, factor);
	V->enforceBoundary();
//	HinaPE::ProjectToNonDivergent(V, &MARKER);
	V->projectToNonDivergent();
	V->enforceBoundary();
	V->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.f);
	V->enforceBoundary();


	// Density Step
	EmitSource(D->getField(), S->getField());
//	HinaPE::Diffuse(D->getField(), &MARKER, factor);
	D->enforceBoundary();
	D->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.f);
	D->enforceBoundary();


	// Temperature Step
	EmitSource(T->getField(), S->getField());
	HinaPE::Diffuse(T->getField(), &MARKER, factor);
	D->enforceBoundary();
	T->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.f);
	T->enforceBoundary();

	return true;
}
