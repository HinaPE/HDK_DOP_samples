#ifndef HDK_DOP_SAMPLES_SIM_DFSPH_PARTICLES_H
#define HDK_DOP_SAMPLES_SIM_DFSPH_PARTICLES_H

#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_Geometry.h>

#define GETSET_FUNCS_BOOL(PRM_NAME) GETSET_DATA_FUNCS_B(#PRM_NAME, PRM_NAME)
#define GETSET_FUNCS_INT(PRM_NAME) GETSET_DATA_FUNCS_I(#PRM_NAME, PRM_NAME)
#define GETSET_FUNCS_FLOAT(PRM_NAME) GETSET_DATA_FUNCS_F(#PRM_NAME, PRM_NAME)
#define GETSET_FUNCS_STRING(PRM_NAME) GETSET_DATA_FUNCS_S(#PRM_NAME, PRM_NAME)
#define GETSET_FUNCS_VECTOR3(PRM_NAME) GETSET_DATA_FUNCS_V3(#PRM_NAME, PRM_NAME)

#define PARAMETER_BOOL(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_TOGGLE, 1, &NAME, &Default##NAME);
#define PARAMETER_INT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_INT, 1, &NAME, &Default##NAME);
#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_STRING(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(0, DEFAULT_VALUE);PRMs.emplace_back(PRM_STRING, 1, &NAME, &Default##NAME);
#define PARAMETER_VECTOR_INT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME);static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__};PRMs.emplace_back(PRM_INT, SIZE, &NAME, Default##NAME.data());
#define PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME);static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__};PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());

#define GEOMETRY_POINT_ATTRIBUTE_INT(ATTRIBUTE_NAME) gdp->addIntTuple(GA_ATTRIB_POINT, ATTRIBUTE_NAME, 1, GA_Defaults(0))->setTypeInfo(GA_TYPE_VOID);
#define GEOMETRY_POINT_ATTRIBUTE_FLOAT(ATTRIBUTE_NAME) gdp->addFloatTuple(GA_ATTRIB_POINT, ATTRIBUTE_NAME, 1, GA_Defaults(0))->setTypeInfo(GA_TYPE_VOID);
#define GEOMETRY_POINT_ATTRIBUTE_STRING(ATTRIBUTE_NAME) gdp->addStringTuple(GA_ATTRIB_POINT, ATTRIBUTE_NAME, 1)->setTypeInfo(GA_TYPE_VOID);
#define GEOMETRY_POINT_ATTRIBUTE_VECTOR_FLOAT_3(ATTRIBUTE_NAME) gdp->addFloatTuple(GA_ATTRIB_POINT, ATTRIBUTE_NAME, 3, GA_Defaults(0))->setTypeInfo(GA_TYPE_VECTOR);

#define SIM_ATTRIBUTE_NAME_V "v"
#define SIM_ATTRIBUTE_NAME_M "m"


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

#endif //HDK_DOP_SAMPLES_SIM_DFSPH_PARTICLES_H
