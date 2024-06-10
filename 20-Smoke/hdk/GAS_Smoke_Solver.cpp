#include "GAS_Smoke_Solver.h"

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
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);

#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_I(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addIntTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleI NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_I(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addIntTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleI NAME##_handle(NAME##_attr);
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

void EmitSourcePartial(SIM_RawField *TARGET, const SIM_RawField *S, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	UT_Interrupt *boss = UTgetInterrupt();
	vit.setArray(TARGET->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3 pos;
		TARGET->indexToPos(vit.x(), vit.y(), vit.z(), pos);

		fpreal value = S->getValue(pos);
		if (value > std::numeric_limits<float>::epsilon())
		{
			value += S->getValue(pos);
			vit.setValue(value);
		}
	}
}
THREADED_METHOD2(, TARGET->shouldMultiThread(), EmitSource, SIM_RawField*, TARGET, const SIM_RawField *, S);

void ComputeBuoyancyPartial(SIM_RawField *V_Y, const SIM_RawField *D, const SIM_RawField *T, const fpreal T_AVG, const fpreal DT, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	UT_Interrupt *boss = UTgetInterrupt();
	vit.setArray(V_Y->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	constexpr fpreal BuoyancyDensityFactor = -0.000625f;
	constexpr fpreal BuoyancyTemperatureFactor = 0.1f;

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
THREADED_METHOD5(, V_Y->shouldMultiThread(), ComputeBuoyancy, SIM_RawField*, V_Y, const SIM_RawField *, D, const SIM_RawField *, T, const fpreal, T_AVG, const fpreal, DT);

void ForwardDiffusionPartial(SIM_RawField *TARGET, const SIM_RawField *S, const fpreal DT, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	UT_Interrupt *boss = UTgetInterrupt();
	vit.setArray(TARGET->fieldNC());
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	constexpr fpreal DiffusionFactor = 0.01f;

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		UT_Vector3 pos;
		TARGET->indexToPos(vit.x(), vit.y(), vit.z(), pos);
//
		fpreal v_s = S->getValue(pos);
		fpreal v_t = v_s + DiffusionFactor * DT * S->getLaplacian(pos);
		vit.setValue(v_t);
	}
}
THREADED_METHOD3(, TARGET->shouldMultiThread(), ForwardDiffusion, SIM_RawField *, TARGET, const SIM_RawField *, S, const fpreal, DT);

bool GAS_Smoke_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_ScalarField *T = getScalarField(obj, GAS_NAME_TEMPERATURE);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_SOURCE);
	SIM_ScalarField *C = getScalarField(obj, GAS_NAME_COLLISION);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);

//	if (!G || !D || !T || !C || !V || !S)
//	{
//		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
//		return false;
//	}

//	SIM_GeometryAutoWriteLock lock(G);
//	GU_Detail &gdp = lock.getGdp();

//	EmitSource(D->getField(), S->getField());
//	EmitSource(T->getField(), S->getField());
//	ComputeBuoyancy(V->getYField(), D->getField(), T->getField(), T->getField()->average(), timestep);
////	SIM_RawField VX_Copy(*V->getXField());
////	SIM_RawField VY_Copy(*V->getYField());
////	SIM_RawField VZ_Copy(*V->getZField());
////	ForwardDiffusion(V->getXField(), &VX_Copy, timestep);
////	ForwardDiffusion(V->getYField(), &VY_Copy, timestep);
////	ForwardDiffusion(V->getZField(), &VZ_Copy, timestep);
//	V->projectToNonDivergent();
//	V->enforceBoundary();
//	V->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
//	V->enforceBoundary();
//	D->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
//	D->enforceBoundary();
//	T->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
//	T->enforceBoundary();

	UT_VoxelArrayF &V_D = *D->getField()->fieldNC();

	V_D.setValue(0, 0, 0, 2.0f);
	printf("Before V_D(0, 0, 0) = %f\n", V_D.getValue(0, 0, 0));
	printf("Before V_D(0, 0, 0) = %f\n", V_D(0, 0, 0));
	V_D.setValue(0, 0, 0, 1.0f);
	printf("After V_D(0, 0, 0) = %f\n", V_D.getValue(0, 0, 0));
	printf("After V_D(0, 0, 0) = %f\n", V_D(0, 0, 0));

	SIM_RawField &R_D = *D->getField();

	return true;
}
