#include "GAS_Smoke_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_FieldUtils.h>
#include <GAS/GAS_ProjectNonDivergent.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);

#define GLOBAL_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);

void GAS_Smoke_Solver::initializeSubclass()
{
	SIM_Data::initializeSubclass();
	this->V_X_tmp = std::make_shared<SIM_RawField>();
	this->V_Y_tmp = std::make_shared<SIM_RawField>();
	this->V_Z_tmp = std::make_shared<SIM_RawField>();
	this->D_tmp = std::make_shared<SIM_RawField>();
	this->T_tmp = std::make_shared<SIM_RawField>();
}
void GAS_Smoke_Solver::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
	this->V_X_tmp = ((GAS_Smoke_Solver *) source)->V_X_tmp;
	this->V_Y_tmp = ((GAS_Smoke_Solver *) source)->V_Y_tmp;
	this->V_Z_tmp = ((GAS_Smoke_Solver *) source)->V_Z_tmp;
	this->D_tmp = ((GAS_Smoke_Solver *) source)->D_tmp;
	this->T_tmp = ((GAS_Smoke_Solver *) source)->T_tmp;
}
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

void advectPartial(SIM_RawField *density, const SIM_RawField *velocity, const int axis, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	UT_Interrupt *boss = UTgetInterrupt();
	vit.setConstArray(velocity->field());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		vit.setValue(1);
	}

	UT_Vector3I backward_face = SIM::FieldUtils::faceToCellMap({vit.x(), vit.y(), vit.z()}, axis, 0);
	UT_Vector3I forward_face = SIM::FieldUtils::faceToCellMap({vit.x(), vit.y(), vit.z()}, axis, 1);
	SIM::FieldUtils::getFieldValue(*velocity, forward_face);
}
THREADED_METHOD3(, true, advect, SIM_RawField *, density, const SIM_RawField *, velocity, const int, axis);

void ComputeBuoyancyPartial(SIM_RawField *V_Y, const SIM_RawField *D, const SIM_RawField *T, const fpreal T_AVG, const fpreal DT, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	UT_Interrupt *boss = UTgetInterrupt();
	vit.setArray(V_Y->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	constexpr fpreal BuoyancyDensityFactor = -0.000625f;
	constexpr fpreal BuoyancyTemperatureFactor = 1.f;

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3 pos;
		V_Y->indexToPos(vit.x(), vit.y(), vit.z(), pos);

		fpreal d = D->getValue(pos);
		fpreal t = T->getValue(pos);
		fpreal v = vit.getValue();
		v += DT * (BuoyancyDensityFactor * d + BuoyancyTemperatureFactor * (t - T_AVG));
		vit.setValue(v);
	}
}
THREADED_METHOD5(, true, ComputeBuoyancy, SIM_RawField*, V_Y, const SIM_RawField *, D, const SIM_RawField *, T, const fpreal, T_AVG, const fpreal, DT);

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

//	SIM_GeometryAutoWriteLock lock(G);
//	GU_Detail &gdp = lock.getGdp();
//
//	UT_Vector3 ScalarFieldOrigin = D->getOrig();
//	UT_Vector3 RawFieldOrigin = D->getField()->getOrig();
//	UT_Vector3 ScalarFieldRes = D->getTotalVoxelRes();
//	UT_Vector3I RawFieldRes = D->getField()->getVoxelRes();
//	UT_Vector3I FieldRes = D->getField()->field()->getVoxelRes();
//	GLOBAL_ATTRIBUTE_V3(ScalarFieldOrigin)
//	ScalarFieldOrigin_handle.set(0, ScalarFieldOrigin);
//	GLOBAL_ATTRIBUTE_V3(RawFieldOrigin)
//	RawFieldOrigin_handle.set(0, RawFieldOrigin);
//	GLOBAL_ATTRIBUTE_V3(ScalarFieldRes)
//	ScalarFieldRes_handle.set(0, ScalarFieldRes);
//	GLOBAL_ATTRIBUTE_V3(RawFieldRes)
//	RawFieldRes_handle.set(0, RawFieldRes);
//	GLOBAL_ATTRIBUTE_V3(FieldRes)
//	FieldRes_handle.set(0, FieldRes);
//
//	UT_Vector3 ScalarField000;
//	UT_Vector3 RawFieldCell000;
//	UT_Vector3 RawField000;
//	UT_Vector3 Field000;
//	D->indexToPos(0, 0, 0, ScalarField000);
//	D->getField()->cellIndexToPos(0, 0, 0, RawFieldCell000);
//	D->getField()->indexToPos(0, 0, 0, RawField000);
//	D->getField()->field()->indexToPos(0, 0, 0, Field000);
//
//	GLOBAL_ATTRIBUTE_V3(ScalarField000)
//	ScalarField000_handle.set(0, ScalarField000);
//	GLOBAL_ATTRIBUTE_V3(RawFieldCell000)
//	RawFieldCell000_handle.set(0, RawFieldCell000);
//	GLOBAL_ATTRIBUTE_V3(RawField000)
//	RawField000_handle.set(0, RawField000);
//	GLOBAL_ATTRIBUTE_V3(Field000)
//	Field000_handle.set(0, Field000);

//	V_X_tmp->match(*V->getXField());
//	V_X_tmp->makeConstant(0.0f);
//	V_Y_tmp->match(*V->getYField());
//	V_Y_tmp->makeConstant(0.0f);
//	V_Z_tmp->match(*V->getZField());
//	V_Z_tmp->makeConstant(0.0f);
//	D_tmp->match(*D->getField());
//	D_tmp->makeConstant(0.0f);
//	T_tmp->match(*T->getField());
//	T_tmp->makeConstant(0.0f);

	ComputeBuoyancy(V->getYField(), D->getField(), T->getField(), T->getField()->average(), timestep);
	V->projectToNonDivergent();
	V->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
	D->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
	T->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
	V->enforceBoundary();
	D->enforceBoundary();
	T->enforceBoundary();

	return true;
}
