#include "GAS_FMS_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_RootData.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GuidePerObject.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>

#include "src/FastMassSpring.h"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define GUIDE_PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME); static PRM_Default Default##NAME(DEFAULT_VALUE); GUIDE_PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);

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
// Return a shared guide so that we only have to build our geometry once. But set the displayonce flag to false so that we can set a different transform for each object.
SIM_Guide *GAS_FMS_Solver::createGuideObjectSubclass() const { return new SIM_GuidePerObject(this); }
void GAS_FMS_Solver::buildGuideGeometrySubclass(const SIM_RootData &root, const SIM_Options &options, const GU_DetailHandle &gdh, UT_DMatrix4 *xform, const SIM_Time &t) const
{
	if (gdh.isNull()) return;
	GU_DetailHandleAutoWriteLock gdl(gdh);
	GU_Detail *gdp = gdl.getGdp();

	if (this->ImplSIMD == nullptr)
		return;

	if (this->_x.length() == 0 || this->_v.length() == 0 || this->_f.length() == 0)
		return;

	auto scale = getScale(options);

	exint size = this->_f.length() / 3;
	for (int i = 0; i < size; ++i)
	{
		// Line
		GA_Offset p1 = gdp->appendPoint();
		GA_Offset p2 = gdp->appendPoint();
		UT_Vector3 pos1 = UT_Vector3(this->_x(3 * i + 0), this->_x(3 * i + 1), this->_x(3 * i + 2));
		UT_Vector3 pos2 = pos1 + scale * UT_Vector3(this->_f(3 * i + 0), this->_f(3 * i + 1), this->_f(3 * i + 2));
		gdp->setPos3(p1, pos1);
		gdp->setPos3(p2, pos2);
		GEO_PrimPoly *line = GEO_PrimPoly::build(gdp, 2, true, false);
		line->setPointOffset(0, p1);
		line->setPointOffset(1, p2);

		// Arrow
		GA_Offset p3 = gdp->appendPoint();
		GA_Offset p4 = gdp->appendPoint();
		GA_Offset p5 = gdp->appendPoint();
		GA_Offset p6 = gdp->appendPoint();
		UT_Vector3 dir = (pos2 - pos1);
		dir.normalize();
		UT_Vector3 tmp = UT_Vector3{1, 0, 0};
		if (tmp == dir)
			tmp = UT_Vector3{0, 1, 0};
		UT_Vector3 right = dir;
		right.cross(tmp);
		right.normalize();
		UT_Vector3 up = right;
		up.cross(dir);
		dir *= 0.01;
		right *= 0.005;
		up *= 0.005;
		UT_Vector3 pos3 = pos2 - dir + right;
		UT_Vector3 pos4 = pos2 - dir - right;
		UT_Vector3 pos5 = pos2 - dir + up;
		UT_Vector3 pos6 = pos2 - dir - up;
		gdp->setPos3(p3, pos3);
		gdp->setPos3(p4, pos4);
		gdp->setPos3(p5, pos5);
		gdp->setPos3(p6, pos6);
		GEO_PrimPoly *arrow1 = GEO_PrimPoly::build(gdp, 2, true, false);
		arrow1->setPointOffset(0, p2);
		arrow1->setPointOffset(1, p3);
		GEO_PrimPoly *arrow2 = GEO_PrimPoly::build(gdp, 2, true, false);
		arrow2->setPointOffset(0, p2);
		arrow2->setPointOffset(1, p4);
		GEO_PrimPoly *arrow3 = GEO_PrimPoly::build(gdp, 2, true, false);
		arrow3->setPointOffset(0, p2);
		arrow3->setPointOffset(1, p5);
		GEO_PrimPoly *arrow4 = GEO_PrimPoly::build(gdp, 2, true, false);
		arrow4->setPointOffset(0, p2);
		arrow4->setPointOffset(1, p6);
	}

	(*xform) = this->_xform;
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

	static std::vector<PRM_Template> GUIDE_PRMs;
	GUIDE_PARAMETER_FLOAT(GuideScale, 1.f)
	GUIDE_PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   DOP_NAME,
								   DOP_ENGLISH,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	DESC.setGuideTemplates(GUIDE_PRMs.data());
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

	// Fetch Force
	GA_RWAttributeRef f_attr = gdp.findPointAttribute(SIM_NAME_FORCE);
	if (f_attr.isInvalid())
		f_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, SIM_NAME_FORCE, 3, GA_Defaults(0));
	GA_RWHandleV3 f_handle(f_attr);

	// Add Gravity
	UT_Vector3F g = {0, 0, 0};
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
		std::vector<std::array<float, 3>> init_x;
		std::set<std::pair<size_t, size_t>> springs;
		GEO_Primitive *prim;
		init_x.resize(gdp.getNumPoints());
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				UT_Vector3F pos = p_handle.get(pt_off);
				init_x[pt_idx] = {pos.x(), pos.y(), pos.z()};
			}
		GA_FOR_ALL_PRIMITIVES(&gdp, prim)
		{
			GEO_PrimPoly *poly = static_cast<GEO_PrimPoly *>(prim);
			if (poly->getVertexCount() == 3)
			{
				GA_Index v1 = poly->getPointIndex(0);
				GA_Index v2 = poly->getPointIndex(1);
				GA_Index v3 = poly->getPointIndex(2);
				v1 < v2 ? springs.insert({v1, v2}) : springs.insert({v2, v1});
				v2 < v3 ? springs.insert({v2, v3}) : springs.insert({v3, v2});
				v3 < v1 ? springs.insert({v3, v1}) : springs.insert({v1, v3});
			} else if (poly->getVertexCount() == 4)
			{
				GA_Index v1 = poly->getPointIndex(0);
				GA_Index v2 = poly->getPointIndex(1);
				GA_Index v3 = poly->getPointIndex(2);
				GA_Index v4 = poly->getPointIndex(3);
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
		ImplSIMD->build(init_x, springs);
	}
	ImplSIMD->Param.gravity[0] = g[0];
	ImplSIMD->Param.gravity[1] = g[1];
	ImplSIMD->Param.gravity[2] = g[2];
	{
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				UT_Vector3F pos = p_handle.get(pt_off);
				ImplSIMD->Cloth->x(3 * pt_idx + 0) = pos.x();
				ImplSIMD->Cloth->x(3 * pt_idx + 1) = pos.y();
				ImplSIMD->Cloth->x(3 * pt_idx + 2) = pos.z();

				UT_Vector3F vel = v_handle.get(pt_off);
				ImplSIMD->Cloth->v(3 * pt_idx + 0) = vel.x();
				ImplSIMD->Cloth->v(3 * pt_idx + 1) = vel.y();
				ImplSIMD->Cloth->v(3 * pt_idx + 2) = vel.z();
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
				p_handle.set(pt_off, UT_Vector3F(ImplSIMD->Cloth->x(3 * pt_idx + 0), ImplSIMD->Cloth->x(3 * pt_idx + 1), ImplSIMD->Cloth->x(3 * pt_idx + 2)));
				v_handle.set(pt_off, UT_Vector3F(ImplSIMD->Cloth->v(3 * pt_idx + 0), ImplSIMD->Cloth->v(3 * pt_idx + 1), ImplSIMD->Cloth->v(3 * pt_idx + 2)));
				f_handle.set(pt_off, UT_Vector3F(ImplSIMD->Cloth->f(3 * pt_idx + 0), ImplSIMD->Cloth->f(3 * pt_idx + 1), ImplSIMD->Cloth->f(3 * pt_idx + 2)));
			}
	}

	exint sz = ImplSIMD->Cloth->x.length();
	this->_x.init(0, sz - 1);
	this->_v.init(0, sz - 1);
	this->_f.init(0, sz - 1);
	this->_x = ImplSIMD->Cloth->x;
	this->_v = ImplSIMD->Cloth->v;
	this->_f = ImplSIMD->Cloth->f;

	UT_DMatrix4 local1, local2;
	G->getTransform(local1);
	const SIM_Position *pos = obj->getPositionForGeometry(SIM_GEOMETRY_DATANAME);
	pos->getTransform(local2);
	_xform = local1 * local2;
	
	return true;
}
