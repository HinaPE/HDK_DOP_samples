#include "GAS_FMS_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>

#include "src/FastMassSpring.h"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define SIM_NAME_ACCELERATION "a"

void GAS_FMS_Solver::initializeSubclass()
{
	SIM_Data::initializeSubclass();
	this->ImplSIMD = nullptr;
}
void GAS_FMS_Solver::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
	this->ImplSIMD = static_cast<const GAS_FMS_Solver *>(source)->ImplSIMD;
}
const SIM_DopDescription *GAS_FMS_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	static std::array<PRM_Name, 8> NumericMethods = {
			PRM_Name("0", "Explicit Euler"),
			PRM_Name("1", "Explicit Symplectic"),
			PRM_Name("2", "Implicit Euler Baraff Witkin"),
			PRM_Name("3", "Integration Gradient Descent"),
			PRM_Name("4", "Integration Newton Descent"),
			PRM_Name("5", "Integration Newton Descent PCG"),
			PRM_Name("6", "Integration Local Global"),
			PRM_Name(nullptr),};
	static PRM_Name NumericMethodsName("NumericMethods", "NumericMethods");
	static PRM_Default NumericMethodsNameDefault(1);
	static PRM_ChoiceList CLNumericMethods(PRM_CHOICELIST_SINGLE, NumericMethods.data());
	PRMs.emplace_back(PRM_ORD, 1, &NumericMethodsName, &NumericMethodsNameDefault, &CLNumericMethods);
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
	GA_Offset pt_off;

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

	// Update Parameters
	if (!ImplSIMD)
	{
		ImplSIMD = std::make_shared<HinaPE::SIMD::FastMassSpring>();
		std::set<std::pair<size_t, size_t>> springs;
		GEO_Primitive *prim;
		GA_FOR_ALL_PRIMITIVES(&gdp, prim)
		{
			GEO_PrimPoly *poly = static_cast<GEO_PrimPoly *>(prim);
			if (poly->getVertexCount() == 3)
			{
				GA_Index v1 = poly->getVertexIndex(0);
				GA_Index v2 = poly->getVertexIndex(1);
				GA_Index v3 = poly->getVertexIndex(2);
				v1 < v2 ? springs.insert({v1, v2}) : springs.insert({v2, v1});
				v2 < v3 ? springs.insert({v2, v3}) : springs.insert({v3, v2});
				v3 < v1 ? springs.insert({v3, v1}) : springs.insert({v1, v3});
			} else if (poly->getVertexCount() == 4)
			{
				GA_Index v1 = poly->getVertexIndex(0);
				GA_Index v2 = poly->getVertexIndex(1);
				GA_Index v3 = poly->getVertexIndex(2);
				GA_Index v4 = poly->getVertexIndex(3);
				v1 < v2 ? springs.insert({v1, v2}) : springs.insert({v2, v1});
				v2 < v3 ? springs.insert({v2, v3}) : springs.insert({v3, v2});
				v3 < v4 ? springs.insert({v3, v4}) : springs.insert({v4, v3});
				v4 < v1 ? springs.insert({v4, v1}) : springs.insert({v1, v4});
				v1 < v3 ? springs.insert({v1, v3}) : springs.insert({v3, v1});
//				v2 < v4 ? springs.insert({v2, v4}) : springs.insert({v4, v2});
			} else
			{
				addError(obj, SIM_MESSAGE, "Invalid Primitive", UT_ERROR_FATAL);
				return false;
			}
		}
		ImplSIMD->build(gdp.getNumPoints(), springs);
	}
	ImplSIMD->Param.gravity[0] = g[0];
	ImplSIMD->Param.gravity[1] = g[1];
	ImplSIMD->Param.gravity[2] = g[2];
	{
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				UT_Vector3 pos = p_handle.get(pt_off);
				ImplSIMD->Cloth->x[pt_idx] = {pos.x(), pos.y(), pos.z()};

				UT_Vector3 vel = v_handle.get(pt_off);
				ImplSIMD->Cloth->v[pt_idx] = {vel.x(), vel.y(), vel.z()};
			}
	}

	switch (getNumericMethods())
	{
		case 0:
			ImplSIMD->solve_explicit_euler(timestep);
			break;
		case 1:
			ImplSIMD->solve_explicit_symplectic(timestep);
			break;
		case 2:
			ImplSIMD->solve_implicit_euler(timestep);
			break;
		case 3:
			ImplSIMD->solve_gradient_descent(timestep);
			break;
		case 4:
			ImplSIMD->solve_newton_descent(timestep);
			break;
		case 5:
			ImplSIMD->solve_newton_descent_pcg(timestep);
			break;
		case 6:
			ImplSIMD->solve_local_global(timestep);
			break;
		default:
		{
			addError(obj, SIM_MESSAGE, "Invalid Numeric Method", UT_ERROR_FATAL);
			return false;
		}
	}

	{
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				p_handle.set(pt_off, UT_Vector3(ImplSIMD->Cloth->x[pt_idx][0], ImplSIMD->Cloth->x[pt_idx][1], ImplSIMD->Cloth->x[pt_idx][2]));
				v_handle.set(pt_off, UT_Vector3(ImplSIMD->Cloth->v[pt_idx][0], ImplSIMD->Cloth->v[pt_idx][1], ImplSIMD->Cloth->v[pt_idx][2]));
				a_handle.set(pt_off, UT_Vector3(ImplSIMD->Cloth->a[pt_idx][0], ImplSIMD->Cloth->a[pt_idx][1], ImplSIMD->Cloth->a[pt_idx][2]));
			}
	}

	return true;
}
