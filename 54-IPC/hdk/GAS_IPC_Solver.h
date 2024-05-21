#ifndef GAS_IPC_SOLVER_H
#define GAS_IPC_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_IPC_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "IPC_Solver";
	inline static const char *DOP_ENGLISH = "IPC Solver";
	inline static const char *DATANAME = "IPC_Solver";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_DATA_FUNCS_I("SubSteps", SubSteps)

protected:
	explicit GAS_IPC_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_IPC_Solver, GAS_SubSolver, "This is a IPC Solver.", getDopDescription());
};

#endif //GAS_IPC_SOLVER_H
