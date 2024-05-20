#include "GAS_Smoke_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <SIM/SIM_ScalarField.h>
#include <GAS/GAS_ProjectNonDivergent.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <SIM/SIM_FieldUtils.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);

#define GLOBAL_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);

const SIM_DopDescription *GAS_Smoke_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	ACTIVATE_GAS_SOURCE
	ACTIVATE_GAS_DENSITY
	ACTIVATE_GAS_TEMPERATURE
	ACTIVATE_GAS_COLLISION
	ACTIVATE_GAS_VELOCITY
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

struct SolverImpl
{
	THREADED_METHOD1(SolverImpl, true, advect, SIM_RawField *, density);
	void advectPartial(SIM_RawField *density, const UT_JobInfo &info)
	{
		UT_VoxelArrayIteratorF vit;
		UT_Interrupt *boss = UTgetInterrupt();
		vit.setArray(density->fieldNC());
		vit.setCompressOnExit(true);
		vit.setPartialRange(info.job(), info.numJobs());

		for (vit.rewind(); !vit.atEnd(); vit.advance())
		{
			vit.setValue(1);
		}
	}
};

bool GAS_Smoke_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_ScalarField *T = getScalarField(obj, GAS_NAME_TEMPERATURE);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_SOURCE);
	SIM_ScalarField *C = getScalarField(obj, GAS_NAME_COLLISION);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);

//	if (!G || !D || !T || !C || !V || !S)
	if (!G || !D || !T)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	SIM_GeometryAutoWriteLock lock(G);
	GU_Detail &gdp = lock.getGdp();

	UT_Vector3 ScalarFieldOrigin = D->getOrig();
	UT_Vector3 RawFieldOrigin = D->getField()->getOrig();
	UT_Vector3 ScalarFieldRes = D->getTotalVoxelRes();
	UT_Vector3I RawFieldRes = D->getField()->getVoxelRes();
	UT_Vector3I FieldRes = D->getField()->field()->getVoxelRes();
	GLOBAL_ATTRIBUTE_V3(ScalarFieldOrigin)
	ScalarFieldOrigin_handle.set(0, ScalarFieldOrigin);
	GLOBAL_ATTRIBUTE_V3(RawFieldOrigin)
	RawFieldOrigin_handle.set(0, RawFieldOrigin);
	GLOBAL_ATTRIBUTE_V3(ScalarFieldRes)
	ScalarFieldRes_handle.set(0, ScalarFieldRes);
	GLOBAL_ATTRIBUTE_V3(RawFieldRes)
	RawFieldRes_handle.set(0, RawFieldRes);
	GLOBAL_ATTRIBUTE_V3(FieldRes)
	FieldRes_handle.set(0, FieldRes);

	UT_Vector3 ScalarField000;
	UT_Vector3 RawFieldCell000;
	UT_Vector3 RawField000;
	UT_Vector3 Field000;
	D->indexToPos(0, 0, 0, ScalarField000);
	D->getField()->cellIndexToPos(0, 0, 0, RawFieldCell000);
	D->getField()->indexToPos(0, 0, 0, RawField000);
	D->getField()->field()->indexToPos(0, 0, 0, Field000);

	GLOBAL_ATTRIBUTE_V3(ScalarField000)
	ScalarField000_handle.set(0, ScalarField000);
	GLOBAL_ATTRIBUTE_V3(RawFieldCell000)
	RawFieldCell000_handle.set(0, RawFieldCell000);
	GLOBAL_ATTRIBUTE_V3(RawField000)
	RawField000_handle.set(0, RawField000);
	GLOBAL_ATTRIBUTE_V3(Field000)
	Field000_handle.set(0, Field000);

	static std::shared_ptr<SolverImpl> Solver = std::make_shared<SolverImpl>();
	Solver->advect(D->getField());

	return true;
}
