#ifndef HDK_DOP_SAMPLES_GAS_DFSPH_SOLVER_H
#define HDK_DOP_SAMPLES_GAS_DFSPH_SOLVER_H

#include <GAS/GAS_SubSolver.h>

#define TARGET_SOLVE_GEOMETRY(PARTICLE_CLASS) static PRM_Name theGeometryName(GAS_NAME_GEOMETRY, PARTICLE_CLASS::DATANAME);static PRM_Default theGeometryNameDefault(0, PARTICLE_CLASS::DATANAME);PRMs.emplace_back(PRM_STRING, 1, &theGeometryName, &theGeometryNameDefault);

class GAS_DFSPH_Solver : public GAS_SubSolver
{
protected:
	explicit GAS_DFSPH_Solver(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_DFSPH_Solver, GAS_SubSolver, "This is a DFSPH Solver.", getDopDescription());

public:
	inline static const bool GEN_NODE = true;
	inline static const char *ENGLISH_NAME = "DFSPH_Solver";
	inline static const char *COMMON_NAME = "DFSPH_Solver";
	inline static const char *DATANAME = "DFSPH_Solver";
	inline static const bool UNIQUE_DATANAME = false;
};

#endif //HDK_DOP_SAMPLES_GAS_DFSPH_SOLVER_H
