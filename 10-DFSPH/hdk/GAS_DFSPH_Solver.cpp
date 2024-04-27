#include "GAS_DFSPH_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GuideShared.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GEO/GEO_PrimPoly.h>

#include "src_tbb/DFSPH.h"
#include "src_simd/DFSPH.h"
#include "src_gpu/DFSPH.cuh"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());
#define GUIDE_PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; GUIDE_PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());

void GAS_DFSPH_Solver::initializeSubclass()
{
	SIM_Data::initializeSubclass();
	this->Impl = nullptr;
}
void GAS_DFSPH_Solver::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
	this->Impl = ((GAS_DFSPH_Solver *) source)->Impl;
}
const SIM_DopDescription *GAS_DFSPH_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	PARAMETER_FLOAT(KernelRadius, 0.04f)
	PARAMETER_VECTOR_FLOAT_N(SolverDomain, 3, 1.f, 1.f, 1.f)
	static std::array<PRM_Name, 4> Backends = {
			PRM_Name("0", "TBB"),
			PRM_Name("1", "SIMD"),
			PRM_Name("2", "CUDA"),
			PRM_Name(nullptr),};
	static PRM_Name BackendsName("Backends", "Backends");
	static PRM_Default BackendsNameDefault(2, "CUDA");
	static PRM_ChoiceList CLBackends(PRM_CHOICELIST_SINGLE, Backends.data());
	PRMs.emplace_back(PRM_ORD, 1, &BackendsName, &BackendsNameDefault, &CLBackends);
	PRMs.emplace_back();

	static std::vector<PRM_Template> GUIDE_PRMs;
	GUIDE_PARAMETER_VECTOR_FLOAT_N(GuideSolverDomain, 3, 1.f, 1.f, 1.f)
	GUIDE_PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   ENGLISH_NAME,
								   COMMON_NAME,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	DESC.setGuideTemplates(GUIDE_PRMs.data());
	setGasDescription(DESC);
	return &DESC;
}
bool GAS_DFSPH_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	if (!G)
		return false;

	SIM_GeometryAutoWriteLock lock(G);
	GU_Detail &gdp = lock.getGdp();

	if (!Impl)
		Impl = std::make_shared<HinaPE::CUDA::DFSPH>(getKernelRadius());
	Impl->resize(gdp.getNumPoints());

	{
		GA_Offset pt_off;
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				UT_Vector3 pos = gdp.getPos3(pt_off);
				Impl->Fluid->x[pt_idx] = {pos.x(), pos.y(), pos.z()};
			}
	}

	Impl->solve(timestep);
	{
		GA_RWAttributeRef v_attr = gdp.findPointAttribute("v");
		if (!v_attr.isValid())
			v_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "v", 3, GA_Defaults(0));
		GA_RWHandleV3 v_handle(v_attr);
		GA_RWAttributeRef rho_attr = gdp.findPointAttribute("rho");
		if (!rho_attr.isValid())
			rho_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "rho", 1, GA_Defaults(0));
		GA_RWHandleF rho_handle(rho_attr);
		GA_RWAttributeRef factor_attr = gdp.findPointAttribute("factor");
		if (!factor_attr.isValid())
			factor_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "factor", 1, GA_Defaults(0));
		GA_RWHandleF factor_handle(factor_attr);
		GA_RWAttributeRef V_attr = gdp.findPointAttribute("V");
		if (!V_attr.isValid())
			V_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "V", 1, GA_Defaults(0));
		GA_RWHandleF V_handle(V_attr);
		GA_RWAttributeRef nn_attr = gdp.findPointAttribute("nn");
		if (!nn_attr.isValid())
			nn_attr = gdp.addIntTuple(GA_ATTRIB_POINT, "nn", 1, GA_Defaults(0));
		GA_RWHandleI nn_handle(nn_attr);
		GA_Offset pt_off;
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				float rho = Impl->Fluid->rho[pt_idx];
				float factor = Impl->Fluid->factor[pt_idx];
				float V = Impl->Fluid->V[pt_idx];
				float nn = Impl->Fluid->nn[pt_idx];
				UT_Vector3 vel = {Impl->Fluid->v[pt_idx].x, Impl->Fluid->v[pt_idx].y, Impl->Fluid->v[pt_idx].z};
				UT_Vector3 pos = {Impl->Fluid->x[pt_idx].x, Impl->Fluid->x[pt_idx].y, Impl->Fluid->x[pt_idx].z};
				rho_handle.set(pt_off, rho);
				factor_handle.set(pt_off, factor);
				V_handle.set(pt_off, V);
				nn_handle.set(pt_off, nn);
				v_handle.set(pt_off, vel);
				gdp.setPos3(pt_off, pos);
			}
	}

	return true;
}
SIM_Guide *GAS_DFSPH_Solver::createGuideObjectSubclass() const { return new SIM_GuideShared(this, false); }
void GAS_DFSPH_Solver::buildGuideGeometrySubclass(const SIM_RootData &root, const SIM_Options &options, const GU_DetailHandle &gdh, UT_DMatrix4 *xform, const SIM_Time &t) const
{
	if (gdh.isNull()) return;
	GU_DetailHandleAutoWriteLock gdl(gdh);
	GU_Detail *gdp = gdl.getGdp();

	UT_Vector3 domain = getSolverDomain(options);

	GA_Offset p1 = gdp->appendPoint();
	GA_Offset p2 = gdp->appendPoint();
	GA_Offset p3 = gdp->appendPoint();
	GA_Offset p4 = gdp->appendPoint();
	GA_Offset p5 = gdp->appendPoint();
	GA_Offset p6 = gdp->appendPoint();
	GA_Offset p7 = gdp->appendPoint();
	GA_Offset p8 = gdp->appendPoint();
	gdp->setPos3(p1, UT_Vector3(-domain.x() / 2, -domain.y() / 2, -domain.z() / 2));
	gdp->setPos3(p2, UT_Vector3(domain.x() / 2, -domain.y() / 2, -domain.z() / 2));
	gdp->setPos3(p3, UT_Vector3(domain.x() / 2, -domain.y() / 2, domain.z() / 2));
	gdp->setPos3(p4, UT_Vector3(-domain.x() / 2, -domain.y() / 2, domain.z() / 2));
	gdp->setPos3(p5, UT_Vector3(-domain.x() / 2, domain.y() / 2, -domain.z() / 2));
	gdp->setPos3(p6, UT_Vector3(domain.x() / 2, domain.y() / 2, -domain.z() / 2));
	gdp->setPos3(p7, UT_Vector3(domain.x() / 2, domain.y() / 2, domain.z() / 2));
	gdp->setPos3(p8, UT_Vector3(-domain.x() / 2, domain.y() / 2, domain.z() / 2));

	GEO_PrimPoly *front = GEO_PrimPoly::build(gdp, 5, true, false);
	GEO_PrimPoly *back = GEO_PrimPoly::build(gdp, 5, true, false);
	GEO_PrimPoly *top = GEO_PrimPoly::build(gdp, 5, true, false);
	GEO_PrimPoly *bottom = GEO_PrimPoly::build(gdp, 5, true, false);
	GEO_PrimPoly *left = GEO_PrimPoly::build(gdp, 5, true, false);
	GEO_PrimPoly *right = GEO_PrimPoly::build(gdp, 5, true, false);
	front->setPointOffset(0, p1);
	front->setPointOffset(1, p2);
	front->setPointOffset(2, p3);
	front->setPointOffset(3, p4);
	front->setPointOffset(4, p1);
	back->setPointOffset(0, p5);
	back->setPointOffset(1, p6);
	back->setPointOffset(2, p7);
	back->setPointOffset(3, p8);
	back->setPointOffset(4, p5);
	top->setPointOffset(0, p1);
	top->setPointOffset(1, p2);
	top->setPointOffset(2, p6);
	top->setPointOffset(3, p5);
	top->setPointOffset(4, p1);
	bottom->setPointOffset(0, p3);
	bottom->setPointOffset(1, p4);
	bottom->setPointOffset(2, p8);
	bottom->setPointOffset(3, p7);
	bottom->setPointOffset(4, p3);
	left->setPointOffset(0, p1);
	left->setPointOffset(1, p4);
	left->setPointOffset(2, p8);
	left->setPointOffset(3, p5);
	left->setPointOffset(4, p1);
	right->setPointOffset(0, p2);
	right->setPointOffset(1, p3);
	right->setPointOffset(2, p7);
	right->setPointOffset(3, p6);
	right->setPointOffset(4, p2);
}
