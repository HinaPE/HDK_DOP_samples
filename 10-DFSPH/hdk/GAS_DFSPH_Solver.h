#ifndef GAS_DFSPH_SOLVER_H
#define GAS_DFSPH_SOLVER_H

#include <GAS/GAS_SubSolver.h>

namespace HinaPE::TBB { struct DFSPH; }
namespace HinaPE::SIMD { struct DFSPH; }
namespace HinaPE::CUDA { struct DFSPH; }

class GAS_DFSPH_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "DFSPH_Solver";
	inline static const char *DOP_ENGLISH = "DFSPH Solver";
	inline static const char *DATANAME = "DFSPH_Solver";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_DATA_FUNCS_F("KernelRadius", KernelRadius)
	GETSET_DATA_FUNCS_I("Backends", Backends)
	GET_GUIDE_FUNC_V3("GuideSolverDomain", SolverDomain, (1, 1, 1));
	std::shared_ptr<HinaPE::TBB::DFSPH> ImplTBB;
	std::shared_ptr<HinaPE::SIMD::DFSPH> ImplSIMD;
	std::shared_ptr<HinaPE::CUDA::DFSPH> ImplCUDA;

protected:
	explicit GAS_DFSPH_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() override;
	void makeEqualSubclass(const SIM_Data *source) override;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	SIM_Guide *createGuideObjectSubclass() const override;
	void buildGuideGeometrySubclass(const SIM_RootData &root, const SIM_Options &options, const GU_DetailHandle &gdh, UT_DMatrix4 *xform, const SIM_Time &t) const override;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_DFSPH_Solver, GAS_SubSolver, "This is a DFSPH Solver.", getDopDescription());
};

#endif //GAS_DFSPH_SOLVER_H
