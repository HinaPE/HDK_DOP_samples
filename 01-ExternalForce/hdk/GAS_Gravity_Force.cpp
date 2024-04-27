#include "GAS_Gravity_Force.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_Position.h>
#include <SIM/SIM_PositionSimple.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());

const SIM_DopDescription *GAS_Gravity_Force::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	PARAMETER_VECTOR_FLOAT_N(Gravity, 3, 0.f, -9.8f, 0.f)
	PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   ENGLISH_NAME,
								   COMMON_NAME,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	setGasDescription(DESC);
	return &DESC;
}
bool GAS_Gravity_Force::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	if (!G)
		return false;

	SIM_GeometryAutoWriteLock lock(G);
	GU_Detail &gdp = lock.getGdp();
	UT_DMatrix4 &trans = G->lockTransform();
	UT_Vector3 pos;
	trans.getTranslates(pos);

	GA_RWAttributeRef v_attr = gdp.findGlobalAttribute(SIM_NAME_VELOCITY);
	if (v_attr.isInvalid())
		v_attr = gdp.addFloatTuple(GA_ATTRIB_GLOBAL, SIM_NAME_VELOCITY, 3, GA_Defaults(0));
	GA_RWHandleV3 v_handle(v_attr);
	UT_Vector3 vel = v_handle.get(0);

	UT_Vector3 a = getGravity();

	vel = vel + timestep * a;
	pos = pos + timestep * vel;

	v_handle.set(0, vel);
	trans.setTranslates(pos);
	G->releaseTransform();

	return true;
}
