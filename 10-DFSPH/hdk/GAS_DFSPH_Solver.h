#ifndef GAS_DFSPH_SOLVER_H
#define GAS_DFSPH_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_DFSPH_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *ENGLISH_NAME = "DFSPH_Solver";
	inline static const char *COMMON_NAME = "DFSPH_Solver";
	inline static const char *DATANAME = "DFSPH_Solver";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_DFSPH_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_DFSPH_Solver, GAS_SubSolver, "This is a DFSPH Solver.", getDopDescription());
};

#endif //GAS_DFSPH_SOLVER_H
