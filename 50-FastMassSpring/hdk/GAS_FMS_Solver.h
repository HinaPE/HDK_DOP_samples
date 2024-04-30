#ifndef GAS_FMS_SOLVER_H
#define GAS_FMS_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_FMS_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "FMS_Solver";
	inline static const char *DOP_ENGLISH = "FMS Solver";
	inline static const char *DATANAME = "FMS_Solver";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_FMS_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_FMS_Solver, GAS_SubSolver, "This is a Fast Mass Spring Solver.", getDopDescription());
};

#endif //GAS_FMS_SOLVER_H
