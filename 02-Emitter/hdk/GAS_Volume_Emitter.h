#ifndef GAS_VOLUME_EMITTER_H
#define GAS_VOLUME_EMITTER_H

#include <GAS/GAS_SubSolver.h>

class GAS_Volume_Emitter : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *ENGLISH_NAME = "Volume_Emitter";
	inline static const char *COMMON_NAME = "Volume_Emitter";
	inline static const char *DATANAME = "Volume_Emitter";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_Volume_Emitter(const SIM_DataFactory *factory) : BaseClass(factory) {}
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_Volume_Emitter, GAS_SubSolver, "This is a Volume Emitter.", getDopDescription());
};

#endif //GAS_VOLUME_EMITTER_H
