#ifndef GAS_GRAVITY_FORCE_H
#define GAS_GRAVITY_FORCE_H

#include <GAS/GAS_SubSolver.h>

class GAS_Gravity_Force : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *ENGLISH_NAME = "Gravity_Force";
	inline static const char *COMMON_NAME = "Gravity_Force";
	inline static const char *DATANAME = "Gravity_Force";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_DATA_FUNCS_V3("Gravity",	Gravity)
	GETSET_DATA_FUNCS_S(GAS_NAME_GEOMETRY, GeometryDATANAME);
	GETSET_DATA_FUNCS_I("TransformTarget", TransformTarget)

protected:
	explicit GAS_Gravity_Force(const SIM_DataFactory *factory) : BaseClass(factory) {}
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_Gravity_Force, GAS_SubSolver, "This is a Gravity Force Solver.", getDopDescription());
};

#endif //GAS_GRAVITY_FORCE_H
