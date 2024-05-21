#ifndef GAS_PBF_SOLVER_H
#define GAS_PBF_SOLVER_H

#include <GAS/GAS_SubSolver.h>

namespace HinaPE::TBB { struct PBF; }

class GAS_PBF_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "PBF_Solver";
	inline static const char *DOP_ENGLISH = "PBF Solver";
	inline static const char *DATANAME = "PBF_Solver";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_DATA_FUNCS_I("Backends", Backends)
	GETSET_DATA_FUNCS_I("SubSteps", SubSteps)
	GETSET_DATA_FUNCS_V3("SolverDomain", Domain)
	GETSET_DATA_FUNCS_F("KernelRadius", KernelRadius)
	GETSET_DATA_FUNCS_F("SurfaceTension", SurfaceTension)
	GETSET_DATA_FUNCS_F("Viscosity", Viscosity)
	GETSET_DATA_FUNCS_B("TopOpen", TopOpen)
	GETSET_DATA_FUNCS_B("DebugInfo", DebugInfo)
	GET_GUIDE_FUNC_V3("GuideSolverDomain", SolverDomain, (1, 1, 1));
	std::shared_ptr<HinaPE::TBB::PBF> ImplTBB;

protected:
	explicit GAS_PBF_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	SIM_Guide *createGuideObjectSubclass() const final;
	void buildGuideGeometrySubclass(const SIM_RootData &root, const SIM_Options &options, const GU_DetailHandle &gdh, UT_DMatrix4 *xform, const SIM_Time &t) const final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_PBF_Solver, GAS_SubSolver, "This is a PBF Solver.", getDopDescription());
};

#endif //GAS_PBF_SOLVER_H
