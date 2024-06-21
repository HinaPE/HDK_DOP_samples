#ifndef GAS_STABLEFLUIDS_H
#define GAS_STABLEFLUIDS_H

#include <GAS/GAS_SubSolver.h>

class GAS_StableFluids : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "StableFluids";
	inline static const char *DOP_ENGLISH = "Stable Fluids";
	inline static const char *DATANAME = "StableFluids";
	inline static const bool UNIQUE_DATANAME = false;

	GETSET_DATA_FUNCS_F("Diffusion", Diffusion)
	GETSET_DATA_FUNCS_F("ImpulseFactor", ImpulseFactor)

protected:
	explicit GAS_StableFluids(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
	DECLARE_STANDARD_GETCASTTOTYPE();
	DECLARE_DATAFACTORY(GAS_StableFluids, GAS_SubSolver, "This is a Stable Fluids Solver.", getDopDescription());
};

#endif //GAS_STABLEFLUIDS_H
