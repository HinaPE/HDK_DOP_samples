#include "GAS_FMS_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define SIM_NAME_ACCELERATION "a"

const SIM_DopDescription *GAS_FMS_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
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
bool GAS_FMS_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	if (!G)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	SIM_GeometryAutoWriteLock lock(G);
	GU_Detail &gdp = lock.getGdp();

	// Fetch Position
	GA_RWHandleV3 p_handle = gdp.getP();

	// Fetch Velocity
	GA_RWAttributeRef v_attr = gdp.findPointAttribute(SIM_NAME_VELOCITY);
	if (v_attr.isInvalid())
		v_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, SIM_NAME_VELOCITY, 3, GA_Defaults(0));
	GA_RWHandleV3 v_handle(v_attr);

	// Fetch Acceleration
	GA_RWAttributeRef a_attr = gdp.findPointAttribute(SIM_NAME_ACCELERATION);
	if (a_attr.isInvalid())
		a_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, SIM_NAME_ACCELERATION, 3, GA_Defaults(0));
	GA_RWHandleV3 a_handle(a_attr);

	// Add Gravity
	UT_Vector3 g = {0, 0, 0};
	SIM_ConstDataArray Forces;
	obj->filterConstSubData(Forces, nullptr, SIM_DataFilterByType("SIM_Force"), SIM_FORCES_DATANAME, SIM_DataFilterNone());
	for (int i = 0; i < Forces.entries(); i++)
	{
		const SIM_ForceGravity *GRV = SIM_DATA_CASTCONST(Forces[i], SIM_ForceGravity);
		if (GRV)
			g += GRV->getGravity();
	};

	GA_Offset pt_off;
	GA_FOR_ALL_PTOFF(&gdp, pt_off)
		{
			UT_Vector3 p = p_handle.get(pt_off);
			UT_Vector3 v = v_handle.get(pt_off);
			UT_Vector3 a = a_handle.get(pt_off);

			a += g;
		}

	return true;
}
