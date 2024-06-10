#include "GAS_KeyFrameSmoke_Solver.h"

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

#include <UT/UT_SparseMatrix.h>

constexpr bool SHOULD_MULTI_THREAD = true;

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_I(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addIntTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleI NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_I(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addIntTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleI NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);

void GAS_KeyFrameSmoke_Solver::initializeSubclass()
{
	SIM_Data::initializeSubclass();
}
void GAS_KeyFrameSmoke_Solver::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
}
const SIM_DopDescription *GAS_KeyFrameSmoke_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	ACTIVATE_GAS_SOURCE
	ACTIVATE_GAS_DENSITY
	ACTIVATE_GAS_TEMPERATURE
	ACTIVATE_GAS_COLLISION
	ACTIVATE_GAS_VELOCITY
	ACTIVATE_GAS_PRESSURE
	ACTIVATE_GAS_DIVERGENCE
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

void EmitSourcePartial(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *SOURCE, float dt, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		fpreal value = vit.getValue();
		value += dt * SOURCE->getValue(vit.x(), vit.y(), vit.z());
		vit.setValue(value);
	}
}
THREADED_METHOD3(, SHOULD_MULTI_THREAD, EmitSource, UT_VoxelArrayF *, TARGET, const UT_VoxelArrayF *, SOURCE, float, dt);

// Bad Diffusion implementation
void BadDiffusePartial(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, float factor, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	const float X0 = ORIGIN->getValue(vit.x(), vit.y(), vit.z());
	for (vit.rewind(); !vit.atEnd(); vit.advance())
		vit.setValue(X0 + factor * (ORIGIN->getValue(vit.x() - 1, vit.y(), vit.z()) +
									ORIGIN->getValue(vit.x() + 1, vit.y(), vit.z()) +
									ORIGIN->getValue(vit.x(), vit.y() - 1, vit.z()) +
									ORIGIN->getValue(vit.x(), vit.y() + 1, vit.z()) +
									ORIGIN->getValue(vit.x(), vit.y(), vit.z() - 1) +
									ORIGIN->getValue(vit.x(), vit.y(), vit.z() + 1) -
									6 * ORIGIN->getValue(vit.x(), vit.y(), vit.z())));
}
THREADED_METHOD3(, SHOULD_MULTI_THREAD, BadDiffuse, UT_VoxelArrayF *, TARGET, const UT_VoxelArrayF *, ORIGIN, float, factor);

// Gauss-Seidel relaxation Diffusion implementation
void DiffuseGaussSeidelPartial(UT_VoxelArrayF *TARGET, const UT_VoxelArrayF *ORIGIN, float factor, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(TARGET);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	const float X0 = ORIGIN->getValue(vit.x(), vit.y(), vit.z());
	for (vit.rewind(); !vit.atEnd(); vit.advance())
		vit.setValue(X0 + factor * (TARGET->getValue(vit.x() - 1, vit.y(), vit.z()) +
									TARGET->getValue(vit.x() + 1, vit.y(), vit.z()) +
									TARGET->getValue(vit.x(), vit.y() - 1, vit.z()) +
									TARGET->getValue(vit.x(), vit.y() + 1, vit.z()) +
									TARGET->getValue(vit.x(), vit.y(), vit.z() - 1) +
									TARGET->getValue(vit.x(), vit.y(), vit.z() + 1) -
									6 * TARGET->getValue(vit.x(), vit.y(), vit.z())));
}
THREADED_METHOD3(, false, DiffuseGaussSeidel, UT_VoxelArrayF *, TARGET, const UT_VoxelArrayF *, ORIGIN, float, factor);

void DivergencePartial(UT_VoxelArrayF *DIVERGENCE, const SIM_VectorField *FLOW, float h, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(DIVERGENCE);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float div = 0;
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		for (int axis: {0, 1, 2})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I face0 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir0);
			UT_Vector3I face1 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir1);
			float v0 = FLOW->getField(axis)->field()->getValue(face0.x(), face0.y(), face0.z());
			float v1 = FLOW->getField(axis)->field()->getValue(face1.x(), face1.y(), face1.z());
			div += (v1 - v0);
		}
		div /= h;
		vit.setValue(div);
	}
}
THREADED_METHOD3(, SHOULD_MULTI_THREAD, Divergence, UT_VoxelArrayF *, DIVERGENCE, const SIM_VectorField *, FLOW, float, h);

void PressureGaussSeidelPartial(UT_VoxelArrayF *PRESSURE, const UT_VoxelArrayF *DIVERGENCE, float h, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(PRESSURE);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float prs = h * h * -DIVERGENCE->getValue(vit.x(), vit.y(), vit.z());
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		for (int axis: {0, 1, 2})
		{
			constexpr int dir0 = 0, dir1 = 1;
			UT_Vector3I face0 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir0);
			UT_Vector3I face1 = SIM::FieldUtils::cellToFaceMap(cell, axis, dir1);
			prs += PRESSURE->getValue(face1.x(), face1.y(), face1.z()) - PRESSURE->getValue(face0.x(), face0.y(), face0.z());
		}
		prs /= 6;
		vit.setValue(prs);
	}
}
THREADED_METHOD3(, SHOULD_MULTI_THREAD, PressureGaussSeidel, UT_VoxelArrayF *, PRESSURE, const UT_VoxelArrayF *, DIVERGENCE, float, h);

void ApplyPressurePartial(UT_VoxelArrayF *FLOW_AXIS, const UT_VoxelArrayF *PRESSURE, float h, int axis, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setArray(FLOW_AXIS);
	vit.setCompressOnExit(true);
	vit.setPartialRange(info.job(), info.numJobs());

	for (vit.rewind(); !vit.atEnd(); vit.advance())
	{
		float v = vit.getValue();
		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		constexpr int dir0 = 0, dir1 = 1;
		UT_Vector3I cell0 = SIM::FieldUtils::faceToCellMap(cell, axis, dir0);
		UT_Vector3I cell1 = SIM::FieldUtils::faceToCellMap(cell, axis, dir1);
		float p0 = PRESSURE->getValue(cell0.x(), cell0.y(), cell0.z());
		float p1 = PRESSURE->getValue(cell1.x(), cell1.y(), cell1.z());
		v -= ((p1 - p0) / h);
		vit.setValue(v);
	}
}
THREADED_METHOD4(, SHOULD_MULTI_THREAD, ApplyPressure, UT_VoxelArrayF *, FLOW_AXIS, const UT_VoxelArrayF *, PRESSURE, float, h, int, axis);

void AdvectPartial(UT_VoxelArrayF *TARGET, const SIM_RawField *ORIGIN, const SIM_VectorField *FLOW, float dt, const UT_JobInfo &info)
{
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(TARGET);
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
THREADED_METHOD4(, SHOULD_MULTI_THREAD, Advect, UT_VoxelArrayF *, TARGET, const SIM_RawField*, SOURCE, const SIM_VectorField *, FLOW, float, dt);

bool GAS_KeyFrameSmoke_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_SOURCE);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);
	SIM_ScalarField *P = getScalarField(obj, GAS_NAME_PRESSURE);
	SIM_ScalarField *Div = getScalarField(obj, GAS_NAME_DIVERGENCE);

	if (!D || !S || !V || !P || !Div)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	float h = D->getField()->getVoxelSize().x();
	constexpr float diff = 0.001f;
	constexpr float visc = 0.001f;
	constexpr size_t gauss_seidel_iterations = 20;

////	 ============================== Velocity Step ==============================
//	// Prepare Buffers
//	UT_VoxelArrayF VX0(*V->getXField()->field());
//	UT_VoxelArrayF VY0(*V->getYField()->field());
//	UT_VoxelArrayF VZ0(*V->getZField()->field());
//	UT_VoxelArrayF VX1(*V->getXField()->field());
//	UT_VoxelArrayF VY1(*V->getYField()->field());
//	UT_VoxelArrayF VZ1(*V->getZField()->field());
//
//	// Apply External Forces
////	EmitSource(&VX0, F->getXField()->field(), timestep);
////	EmitSource(&VY0, F->getYField()->field(), timestep);
////	EmitSource(&VZ0, F->getZField()->field(), timestep);
//
//	// Diffuse Velocity
//	for (size_t i = 0; i < gauss_seidel_iterations; i++)
//		DiffuseGaussSeidelNoThread(&VX1, &VX0, timestep * visc / (h * h));
//	for (size_t i = 0; i < gauss_seidel_iterations; i++)
//		DiffuseGaussSeidelNoThread(&VY1, &VY0, timestep * visc / (h * h));
//	for (size_t i = 0; i < gauss_seidel_iterations; i++)
//		DiffuseGaussSeidelNoThread(&VZ1, &VZ0, timestep * visc / (h * h));
//	V->getXField()->fieldNC()->copyData(VX1);
//	V->getYField()->fieldNC()->copyData(VY1);
//	V->getZField()->fieldNC()->copyData(VZ1);
//
//	// Project Velocity
//	Divergence(Div->getField()->fieldNC(), V, h);
//	P->getField()->fieldNC()->constant(0);
//	for (size_t i = 0; i < gauss_seidel_iterations; i++)
//		PressureGaussSeidelNoThread(P->getField()->fieldNC(), Div->getField()->field(), h);
//	ApplyPressure(V->getXField()->fieldNC(), P->getField()->field(), h, 0);
//	ApplyPressure(V->getYField()->fieldNC(), P->getField()->field(), h, 1);
//	ApplyPressure(V->getZField()->fieldNC(), P->getField()->field(), h, 2);
//
//	// Advect Velocity
//	Advect(&VX0, V->getXField(), V, timestep);
//	Advect(&VY0, V->getYField(), V, timestep);
//	Advect(&VZ0, V->getZField(), V, timestep);
//	V->getXField()->fieldNC()->copyData(VX0);
//	V->getYField()->fieldNC()->copyData(VY0);
//	V->getZField()->fieldNC()->copyData(VZ0);
////	 ============================== Velocity Step ==============================



//	 ============================== Density Step ==============================
	EmitSource(D->getField()->fieldNC(), S->getField()->field(), timestep);
	UT_VoxelArrayF D0(*D->getField()->field());
	UT_VoxelArrayF D1(*D->getField()->field());
//	for (size_t i = 0; i < gauss_seidel_iterations; i++)
		BadDiffuseNoThread(&D1, &D0, timestep * diff / (h * h));
	D->getField()->fieldNC()->copyData(D1);
//	Advect(&D0, D->getField(), V, timestep);
//	D->getField()->fieldNC()->copyData(D0);
//	 ============================== Density Step ==============================

	return true;
}
