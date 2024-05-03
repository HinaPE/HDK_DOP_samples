#ifndef GAS_FMS_SOLVER_H
#define GAS_FMS_SOLVER_H

#include <GAS/GAS_SubSolver.h>

namespace HinaPE::SIMD { struct FastMassSpring; }

class GAS_FMS_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "FMS_Solver";
	inline static const char *DOP_ENGLISH = "FMS Solver";
	inline static const char *DATANAME = "FMS_Solver";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_DATA_FUNCS_I("NumericMethods", NumericMethods)
	GET_GUIDE_FUNC_F("GuideScale", Scale, 1.f)
	GET_GUIDE_FUNC_V3("GuideSolverDomain", SolverDomain, (1, 1, 1));
	std::shared_ptr<HinaPE::SIMD::FastMassSpring> ImplSIMD;

	// ONLY FOR GUIDE GEOMETRY
	mutable UT_VectorF _x;
	mutable UT_VectorF _v;
	mutable UT_VectorF _f;
	mutable UT_DMatrix4 _xform;

protected:
	explicit GAS_FMS_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	SIM_Guide *createGuideObjectSubclass() const override;
	void buildGuideGeometrySubclass(const SIM_RootData &root, const SIM_Options &options, const GU_DetailHandle &gdh, UT_DMatrix4 *xform, const SIM_Time &t) const override;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_FMS_Solver, GAS_SubSolver, "This is a Fast Mass Spring Solver.", getDopDescription());
};

#endif //GAS_FMS_SOLVER_H
