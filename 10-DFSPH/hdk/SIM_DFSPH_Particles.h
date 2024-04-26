#ifndef SIM_DFSPH_PARTICLES_H
#define SIM_DFSPH_PARTICLES_H

#include <SIM/SIM_GeometryCopy.h>
#include "utils.h"

class SIM_DFSPH_Particles : public SIM_GeometryCopy
{
public:
	inline static const bool GEN_NODE = true;
	inline static const char *ENGLISH_NAME = "DFSPH_Particles";
	inline static const char *COMMON_NAME = "DFSPH_Particles";
	inline static const char *DATANAME = "DFSPH_Particles";
	inline static const bool UNIQUE_DATANAME = false;

public:
	GETSET_FUNCS_BOOL(BBB)
	GETSET_FUNCS_INT(III)
	GETSET_FUNCS_FLOAT(FFF)
	GETSET_FUNCS_STRING(STR)
	GETSET_FUNCS_VECTOR3(VEC3I)
	GETSET_FUNCS_VECTOR3(VEC3F)

protected:
	explicit SIM_DFSPH_Particles(const SIM_DataFactory *factory) : BaseClass(factory) {}
	void initializeSubclass() final;
	void makeEqualSubclass(const SIM_Data *source) final;
	GU_ConstDetailHandle getGeometrySubclass() const final;
	mutable GU_DetailHandle _handle;
	static const SIM_DopDescription *getDopDescription();
DECLARE_STANDARD_GETCASTTOTYPE();
DECLARE_DATAFACTORY(SIM_DFSPH_Particles, SIM_GeometryCopy, "This is DFSPH Particles", getDopDescription());
};

#endif //SIM_DFSPH_PARTICLES_H
