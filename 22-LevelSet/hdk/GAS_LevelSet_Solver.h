#ifndef GAS_LEVELSET_SOLVER_H
#define GAS_LEVELSET_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_LevelSet_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "LevelSetSolver";
	inline static const char *DOP_ENGLISH = "LevelSet Solver";
	inline static const char *DATANAME = "LevelSetSolver";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_LevelSet_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
	DECLARE_STANDARD_GETCASTTOTYPE();
	DECLARE_DATAFACTORY(GAS_LevelSet_Solver, GAS_SubSolver, "This is a LevelSet Solver.", getDopDescription());
};

#endif //GAS_LEVELSET_SOLVER_H
