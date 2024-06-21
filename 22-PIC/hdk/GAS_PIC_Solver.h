#ifndef GAS_PIC_SOLVER_H
#define GAS_PIC_SOLVER_H

#include <GAS/GAS_SubSolver.h>

class GAS_PIC_Solver : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "PICSolver";
	inline static const char *DOP_ENGLISH = "PIC Solver";
	inline static const char *DATANAME = "PICSolver";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_PIC_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
	DECLARE_STANDARD_GETCASTTOTYPE();
	DECLARE_DATAFACTORY(GAS_PIC_Solver, GAS_SubSolver, "This is a Particles In Cells Solver.", getDopDescription());
};

#endif //GAS_PIC_SOLVER_H
