#ifndef GAS_RIGID_APPLY_IMPULSE_H
#define GAS_RIGID_APPLY_IMPULSE_H

#include <GAS/GAS_SubSolver.h>

class GAS_Rigid_Apply_Impulse : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *ENGLISH_NAME = "Rigid_Apply_Impulse";
	inline static const char *COMMON_NAME = "Rigid_Apply_Impulse";
	inline static const char *DATANAME = "Rigid_Apply_Impulse";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_DATA_FUNCS_F("Force", Force)
	GETSET_DATA_FUNCS_V3("Position", Position)
	GETSET_DATA_FUNCS_V3("Direction", Direction)

protected:
	explicit GAS_Rigid_Apply_Impulse(const SIM_DataFactory *factory) : BaseClass(factory) {}
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_Rigid_Apply_Impulse, GAS_SubSolver, "This is a Rigid Apply Impulse.", getDopDescription());
};

#endif //GAS_RIGID_APPLY_IMPULSE_H
