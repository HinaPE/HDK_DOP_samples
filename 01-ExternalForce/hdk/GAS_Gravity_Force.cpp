#include "GAS_Gravity_Force.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_Position.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#define ACTIVATE_GAS_GEOMETRY        static PRM_Name GeometryName(GAS_NAME_GEOMETRY, "Geometry"); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);

const SIM_DopDescription *GAS_Gravity_Force::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	PARAMETER_FLOAT(Gravity, -9.8)
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

	const SIM_Position *P = obj->getPositionForGeometry(GAS_NAME_GEOMETRY);
	G->getPositionPath();
	if (P)
	{

	} else
	{
		SIM_GeometryAutoWriteLock lock(G);
		GU_Detail &gdp = lock.getGdp();

		GA_RWAttributeRef v_attr = gdp.findPointAttribute("v");
		if (v_attr.isInvalid())
			v_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "v", 3, GA_Defaults(0));
		GA_RWHandleV3 v_handle(v_attr);

		GA_Offset pt_off;
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				auto pos = gdp.getPos3(pt_off);
				auto vel = v_handle.get(pt_off);
				auto g = getGravity();
				vel.y() += timestep * g;
				pos += timestep * vel;
				v_handle.set(pt_off, vel);
				gdp.setPos3(pt_off, pos);
			}
	}

	return true;
}
