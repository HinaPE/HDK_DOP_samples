#ifndef GAS_KEY_FRAME_SMOKE_SOLVER_H
#define GAS_KEY_FRAME_SMOKE_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_KeyFrameSmoke_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "KeyFrameSmoke_Solver";
	inline static const char *DOP_ENGLISH = "KeyFrameSmoke Solver";
	inline static const char *DATANAME = "KeyFrameSmoke_Solver";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_KeyFrameSmoke_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_KeyFrameSmoke_Solver, GAS_SubSolver, "This is a KeyFrame Smoke Solver.", getDopDescription());
};

#endif //GAS_KEY_FRAME_SMOKE_SOLVER_H
