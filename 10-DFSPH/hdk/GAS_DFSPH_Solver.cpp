#include "GAS_DFSPH_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GuideShared.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GEO/GEO_PrimPoly.h>
#include <GU/GU_NeighbourList.h>

#include "src_tbb/DFSPH.h"
#include "src_simd/DFSPH.h"
#include "src_gpu/DFSPH.cuh"

#include <chrono>

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define PARAMETER_BOOL(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_TOGGLE, 1, &NAME, &Default##NAME);
#define PARAMETER_INT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_INT, 1, &NAME, &Default##NAME);
#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());
#define GUIDE_PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; GUIDE_PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());

void GAS_DFSPH_Solver::initializeSubclass()
{
	SIM_Data::initializeSubclass();
	this->ImplTBB = nullptr;
	this->ImplSIMD = nullptr;
	this->ImplCUDA = nullptr;
}
void GAS_DFSPH_Solver::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
	this->ImplTBB = ((GAS_DFSPH_Solver *) source)->ImplTBB;
	this->ImplSIMD = ((GAS_DFSPH_Solver *) source)->ImplSIMD;
	this->ImplCUDA = ((GAS_DFSPH_Solver *) source)->ImplCUDA;
}
const SIM_DopDescription *GAS_DFSPH_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	static std::array<PRM_Name, 4> Backends = {
			PRM_Name("0", "TBB"),
			PRM_Name("1", "SIMD"),
			PRM_Name("2", "CUDA"),
			PRM_Name(nullptr),};
	static PRM_Name BackendsName("Backends", "Backends");
	static PRM_Default BackendsNameDefault(0);
	static PRM_ChoiceList CLBackends(PRM_CHOICELIST_SINGLE, Backends.data());
	PRMs.emplace_back(PRM_ORD, 1, &BackendsName, &BackendsNameDefault, &CLBackends);
	PARAMETER_INT(SubSteps, 1)
	PARAMETER_VECTOR_FLOAT_N(SolverDomain, 3, 1.f, 1.f, 1.f)
	PARAMETER_FLOAT(KernelRadius, 0.04f)
	PARAMETER_FLOAT(SurfaceTension, 0.01f)
	PARAMETER_FLOAT(Viscosity, 0.01f)
	PARAMETER_BOOL(TopOpen, true)
	PARAMETER_BOOL(DebugInfo, false)
	PRMs.emplace_back();

	static std::vector<PRM_Template> GUIDE_PRMs;
	GUIDE_PARAMETER_VECTOR_FLOAT_N(GuideSolverDomain, 3, 1.f, 1.f, 1.f)
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
bool GAS_DFSPH_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
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
	{
		GA_FOR_ALL_PTOFF(&gdp, pt_off)
			{
				// CHECK POINT INDEX
				GA_Index pt_idx = gdp.pointIndex(pt_off);
				if (pt_idx != pt_off)
					std::cout << "pt_idx != pt_off, pt_idx: " << pt_idx << ", pt_off: " << pt_off << std::endl;
			}
	}
	if (gdp.getNumPoints() == 0)
		return true;

	GA_RWAttributeRef total_time_attr = gdp.findGlobalAttribute("total_time");
	if (!total_time_attr.isValid())
		total_time_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, "total_time", 1, GA_Defaults(0));
	GA_RWHandleF total_time_handle(total_time_attr);
	GA_RWAttributeRef solving_time_attr = gdp.findGlobalAttribute("solving_time");
	if (!solving_time_attr.isValid())
		solving_time_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, "solving_time", 1, GA_Defaults(0));
	GA_RWHandleF solving_time_handle(solving_time_attr);
	GA_RWAttributeRef dt_attr = gdp.findGlobalAttribute("dt");
	if (!dt_attr.isValid())
		dt_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, "dt", 1, GA_Defaults(0));
	GA_RWHandleF dt_handle(dt_attr);
	GA_RWAttributeRef div_iter_attr = gdp.findGlobalAttribute("div_iter");
	if (!div_iter_attr.isValid())
		div_iter_attr = gdp.addIntTuple(GA_ATTRIB_DETAIL, "div_iter", 1, GA_Defaults(0));
	GA_RWHandleI div_iter_handle(div_iter_attr);
	GA_RWAttributeRef prs_iter_attr = gdp.findGlobalAttribute("prs_iter");
	if (!prs_iter_attr.isValid())
		prs_iter_attr = gdp.addIntTuple(GA_ATTRIB_DETAIL, "prs_iter", 1, GA_Defaults(0));
	GA_RWHandleI prs_iter_handle(prs_iter_attr);

	auto start = std::chrono::high_resolution_clock::now();
	for (auto _ = 0; _ < getSubSteps(); ++_)
	{
		switch (getBackends())
		{
			case 0:
			{
				if (!ImplTBB)
					ImplTBB = std::make_shared<HinaPE::TBB::DFSPH>();
				ImplTBB->KERNEL_RADIUS = (float) getKernelRadius();
				ImplTBB->SURFACE_TENSION = (float) getSurfaceTension();
				ImplTBB->VISCOSITY = (float) getViscosity();
				ImplTBB->MaxBound = getDomainF() / 2.f;
				ImplTBB->TOP_OPEN = getTopOpen();
				ImplTBB->solve(timestep / getSubSteps(), &gdp);
				{
					GA_FOR_ALL_PTOFF(&gdp, pt_off)
						{
							GA_Index pt_idx = gdp.pointIndex(pt_off);
							UT_Vector3 pos = {ImplTBB->Fluid->x[pt_idx][0], ImplTBB->Fluid->x[pt_idx][1], ImplTBB->Fluid->x[pt_idx][2]};
							gdp.setPos3(pt_off, pos);
						}

					div_iter_handle.set(0, ImplTBB->DIVERGENCE_ITERS);
					prs_iter_handle.set(0, ImplTBB->PRESSURE_ITERS);
				}
				if (getDebugInfo())
				{
					GA_FOR_ALL_PTOFF(&gdp, pt_off)
						{
							GA_RWAttributeRef v_attr = gdp.findPointAttribute("v");
							if (!v_attr.isValid()) v_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "v", 3, GA_Defaults(0));
							GA_RWHandleV3 v_handle(v_attr);
							GA_RWAttributeRef a_attr = gdp.findPointAttribute("a");
							if (!a_attr.isValid()) a_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "a", 3, GA_Defaults(0));
							GA_RWHandleV3 a_handle(a_attr);
							GA_RWAttributeRef V_attr = gdp.findPointAttribute("V");
							if (!V_attr.isValid()) V_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "V", 1, GA_Defaults(0));
							GA_RWHandleF V_handle(V_attr);
							GA_RWAttributeRef rho_attr = gdp.findPointAttribute("rho");
							if (!rho_attr.isValid()) rho_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "rho", 1, GA_Defaults(0));
							GA_RWHandleF rho_handle(rho_attr);
							GA_RWAttributeRef factor_attr = gdp.findPointAttribute("factor");
							if (!factor_attr.isValid()) factor_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "factor", 1, GA_Defaults(0));
							GA_RWHandleF factor_handle(factor_attr);
							GA_RWAttributeRef nn_attr = gdp.findPointAttribute("nn");
							if (!nn_attr.isValid()) nn_attr = gdp.addIntTuple(GA_ATTRIB_POINT, "nn", 1, GA_Defaults(0));
							GA_RWHandleI nn_handle(nn_attr);

							GA_Index pt_idx = gdp.pointIndex(pt_off);
							UT_Vector3 vel = {ImplTBB->Fluid->v[pt_idx][0], ImplTBB->Fluid->v[pt_idx][1], ImplTBB->Fluid->v[pt_idx][2]};
							UT_Vector3 a = {ImplTBB->Fluid->a[pt_idx][0], ImplTBB->Fluid->a[pt_idx][1], ImplTBB->Fluid->a[pt_idx][2]};
							float V = ImplTBB->Fluid->V[pt_idx];
							float rho = ImplTBB->Fluid->rho[pt_idx];
							float factor = ImplTBB->Fluid->factor[pt_idx];
							float nn = ImplTBB->Fluid->nn[pt_idx];
							v_handle.set(pt_off, vel);
							a_handle.set(pt_off, a);
							V_handle.set(pt_off, V);
							rho_handle.set(pt_off, rho);
							factor_handle.set(pt_off, factor);
							nn_handle.set(pt_off, nn);
						}
				}
			}
				break;
			case 1:
			{
				if (!ImplSIMD)
					ImplSIMD = std::make_shared<HinaPE::SIMD::DFSPH>();
				ImplSIMD->KERNEL_RADIUS = (float) getKernelRadius();
				ImplSIMD->SURFACE_TENSION = (float) getSurfaceTension();
				ImplSIMD->VISCOSITY = (float) getViscosity();
				ImplSIMD->MaxBound = getDomainF() / 2.f;
				ImplSIMD->TOP_OPEN = getTopOpen();
				ImplSIMD->solve(timestep / getSubSteps(), &gdp);
				{
					GA_FOR_ALL_PTOFF(&gdp, pt_off)
						{
							GA_Index pt_idx = gdp.pointIndex(pt_off);
							UT_Vector3 pos = {ImplSIMD->Fluid->x[3 * pt_idx + 0], ImplSIMD->Fluid->x[3 * pt_idx + 1], ImplSIMD->Fluid->x[3 * pt_idx + 2]};
							gdp.setPos3(pt_off, pos);
						}

					div_iter_handle.set(0, ImplSIMD->DIVERGENCE_ITERS);
					prs_iter_handle.set(0, ImplSIMD->PRESSURE_ITERS);
				}
				if (getDebugInfo())
				{
					GA_FOR_ALL_PTOFF(&gdp, pt_off)
						{
							GA_RWAttributeRef v_attr = gdp.findPointAttribute("v");
							if (!v_attr.isValid()) v_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "v", 3, GA_Defaults(0));
							GA_RWHandleV3 v_handle(v_attr);
							GA_RWAttributeRef a_attr = gdp.findPointAttribute("a");
							if (!a_attr.isValid()) a_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "a", 3, GA_Defaults(0));
							GA_RWHandleV3 a_handle(a_attr);
							GA_RWAttributeRef V_attr = gdp.findPointAttribute("V");
							if (!V_attr.isValid()) V_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "V", 1, GA_Defaults(0));
							GA_RWHandleF V_handle(V_attr);
							GA_RWAttributeRef rho_attr = gdp.findPointAttribute("rho");
							if (!rho_attr.isValid()) rho_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "rho", 1, GA_Defaults(0));
							GA_RWHandleF rho_handle(rho_attr);
							GA_RWAttributeRef factor_attr = gdp.findPointAttribute("factor");
							if (!factor_attr.isValid()) factor_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, "factor", 1, GA_Defaults(0));
							GA_RWHandleF factor_handle(factor_attr);
							GA_RWAttributeRef nn_attr = gdp.findPointAttribute("nn");
							if (!nn_attr.isValid()) nn_attr = gdp.addIntTuple(GA_ATTRIB_POINT, "nn", 1, GA_Defaults(0));
							GA_RWHandleI nn_handle(nn_attr);

							GA_Index pt_idx = gdp.pointIndex(pt_off);
							UT_Vector3 vel = {ImplSIMD->Fluid->v[3 * pt_idx + 0], ImplSIMD->Fluid->v[3 * pt_idx + 1], ImplSIMD->Fluid->v[3 * pt_idx + 2]};
							UT_Vector3 a = {ImplSIMD->Fluid->a[3 * pt_idx + 0], ImplSIMD->Fluid->a[3 * pt_idx + 1], ImplSIMD->Fluid->a[3 * pt_idx + 2]};
							float V = ImplSIMD->Fluid->V[pt_idx];
							float rho = ImplSIMD->Fluid->rho[pt_idx];
							float factor = ImplSIMD->Fluid->factor[pt_idx];
							float nn = ImplSIMD->Fluid->nn[pt_idx];
							v_handle.set(pt_off, vel);
							a_handle.set(pt_off, a);
							V_handle.set(pt_off, V);
							rho_handle.set(pt_off, rho);
							factor_handle.set(pt_off, factor);
							nn_handle.set(pt_off, nn);
						}
				}
			}
				break;
			case 2:
			{
				if (!ImplCUDA)
					ImplCUDA = std::make_shared<HinaPE::CUDA::DFSPH>((float) getKernelRadius());
				ImplCUDA->resize(gdp.getNumPoints());

				{
					GA_FOR_ALL_PTOFF(&gdp, pt_off)
						{
							GA_Index pt_idx = gdp.pointIndex(pt_off);
							UT_Vector3 pos = gdp.getPos3(pt_off);
							ImplCUDA->Fluid->x[pt_idx] = {pos.x(), pos.y(), pos.z()};
						}
				}
				ImplCUDA->solve(timestep);
//			{
//				GA_FOR_ALL_PTOFF(&gdp, pt_off)
//					{
//						GA_Index pt_idx = gdp.pointIndex(pt_off);
//
//						UT_Vector3 pos = {ImplCUDA->Fluid->x[pt_idx].x, ImplCUDA->Fluid->x[pt_idx].y, ImplCUDA->Fluid->x[pt_idx].z};
//						UT_Vector3 vel = {ImplCUDA->Fluid->v[pt_idx].x, ImplCUDA->Fluid->v[pt_idx].y, ImplCUDA->Fluid->v[pt_idx].z};
//						UT_Vector3 a = {ImplCUDA->Fluid->a[pt_idx].x, ImplCUDA->Fluid->a[pt_idx].y, ImplCUDA->Fluid->a[pt_idx].z};
//						float rho = ImplCUDA->Fluid->rho[pt_idx];
//						float factor = ImplCUDA->Fluid->factor[pt_idx];
//						float nn = ImplCUDA->Fluid->nn[pt_idx];
//						float V = ImplCUDA->Fluid->V[pt_idx];
//
//						gdp.setPos3(pt_off, pos);
//						v_handle.set(pt_off, vel);
//						a_handle.set(pt_off, a);
//						V_handle.set(pt_off, V);
//						rho_handle.set(pt_off, rho);
//						factor_handle.set(pt_off, factor);
//						nn_handle.set(pt_off, nn);
//					}
//			}
			}
				break;
			default:
			{
				addError(obj, SIM_MESSAGE, "Invalid Backends", UT_ERROR_FATAL);
				return false;
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> duration = end - start;
	solving_time_handle.set(0, (float) duration.count());
	total_time_handle.set(0, total_time_handle.get(0) + (float) duration.count());
	dt_handle.set(0, (float) timestep / getSubSteps());

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
