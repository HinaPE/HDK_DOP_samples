#ifndef GAS_SMOKE_SOLVER_H
#define GAS_SMOKE_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_Smoke_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "Smoke_Solver";
	inline static const char *DOP_ENGLISH = "Smoke Solver";
	inline static const char *DATANAME = "Smoke_Solver";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_Smoke_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_Smoke_Solver, GAS_SubSolver, "This is a Smoke Solver.", getDopDescription());
};

#endif //GAS_SMOKE_SOLVER_H
