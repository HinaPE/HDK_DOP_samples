#ifndef GAS_VOLUME_RENDERER_H
#define GAS_VOLUME_RENDERER_H

#include <GAS/GAS_SubSolver.h>

class GAS_Volume_Renderer : public GAS_SubSolver
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *DOP_NAME = "Volume_Renderer";
	inline static const char *DOP_ENGLISH = "Volume Renderer";
	inline static const char *DATANAME = "Volume_Renderer";
	inline static const bool UNIQUE_DATANAME = false;

protected:
	explicit GAS_Volume_Renderer(const SIM_DataFactory *factory) : BaseClass(factory) {}
	bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep) final;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(GAS_Volume_Renderer, GAS_SubSolver, "This is a Volume Renderer.", getDopDescription());
};

#endif //GAS_VOLUME_RENDERER_H
