#include "GAS_PIC_Solver.h"

#include <SIM/SIM_Engine.h>
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

#include <GAS/GAS_Diffuse.h>
#include <GAS/GAS_Rest.h>

#include <GU/GU_NeighbourList.h>
#include <GAS/GAS_ParticleToSDF.h>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());

inline void parallel_for(size_t n, const std::function<void(size_t)> &f)
{
	UTparallelForEachNumber((int) n, [&](const UT_BlockedRange<int> &range)
	{
		for (size_t i = range.begin(); i != range.end(); ++i) { f(i); }
	});
}

void GAS_PIC_Solver::initializeSubclass() { SIM_Data::initializeSubclass(); }
void GAS_PIC_Solver::makeEqualSubclass(const SIM_Data *source) { SIM_Data::makeEqualSubclass(source); }
const SIM_DopDescription *GAS_PIC_Solver::getDopDescription()
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
	PARAMETER_VECTOR_FLOAT_N(TEST, 3, 0, 0, 0)
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

void Contribute(SIM_RawField *Target, const UT_Vector3 &InPos, const float InValue)
{
	// clamp position to edge
	UT_Vector3 pos = InPos;
	{
		UT_Vector3 low, high;
		Target->indexToPos(0, 0, 0, low);
		Target->indexToPos(Target->getXRes() - 1, Target->getYRes() - 1, Target->getZRes() - 1, high);

		for (size_t axis: {0, 1, 2})
		{
			if (pos[axis] < low[axis])
				pos[axis] = low[axis];
			if (pos[axis] > high[axis])
				pos[axis] = high[axis];
		}
	}

	// compute delta
	int i, j, k;
	float voxel_volume = Target->getVoxelVolume();
	UT_Vector3 delta;
	{
		UT_Vector3 cell_pos;
		Target->posToIndex(pos, i, j, k);
		Target->indexToPos(i, j, k, cell_pos);
		delta = pos - cell_pos;
	}

	// compute weights and transfer value
	std::array<std::array<std::array<float, 2>, 2>, 2> weights;
	for (bool _x: {true, false})
	{
		for (bool _y: {true, false})
		{
			for (bool _z: {true, false})
			{
				float weight = (_x ? std::abs(delta.x()) : Target->getVoxelSize().x() - std::abs(delta.x())) *
				               (_y ? std::abs(delta.y()) : Target->getVoxelSize().y() - std::abs(delta.y())) *
				               (_z ? std::abs(delta.z()) : Target->getVoxelSize().z() - std::abs(delta.z())) /
				               voxel_volume;
				weights[_x][_y][_z] = weight;

				int _i = delta.x() > 0 ? i + _x : i - _x;
				int _j = delta.y() > 0 ? j + _y : j - _y;
				int _k = delta.z() > 0 ? k + _z : k - _z;
				float old_value = Target->field()->getValue(_i, _j, _k);
				Target->fieldNC()->setValue(_i, _j, _k, old_value + InValue * weights[_x][_y][_z]);
			}
		}
	}
}

void P2G(SIM_VectorField *FLOW, const GU_Detail &gdp)
{
	GA_Offset ptoff;
	GA_FOR_ALL_PTOFF(&gdp, ptoff)
	{
	}
}

void G2P(GU_Detail &gdp, const SIM_VectorField *FLOW)
{
	GA_ROHandleV3 p = gdp.getP();
	POINT_ATTRIBUTE_V3(v);

	GA_Offset ptoff;
	GA_FOR_ALL_PTOFF(&gdp, ptoff)
	{
		v_handle.set(ptoff, FLOW->getValue(p.get(ptoff)));
	}
}

bool GAS_PIC_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);

	if (!G || !D || !V)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	SIM_GeometryAutoWriteLock lock(G);
	GU_Detail &gdp = lock.getGdp();
	if (gdp.getNumPoints() == 0)
		return true;

	POINT_ATTRIBUTE_V3(v);
	POINT_ATTRIBUTE_V3(a);

	Contribute(D->getField(), getTEST(), 1);

	return true;
}
