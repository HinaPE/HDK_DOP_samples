#ifndef DFSPH_UTILS_H
#define DFSPH_UTILS_H

#include <GU/GU_Detail.h>

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

inline void RESIZE(GU_Detail &gdp, size_t size)
{
	if (gdp.getNumPoints() > size)
		gdp.deletePoints(GA_Range(gdp.getPointMap(), size, gdp.getNumPoints()));
	if (gdp.getNumPoints() <= size)
		gdp.appendPointBlock(size - gdp.getNumPoints());
}
inline void SYNC_B(GU_Detail &gdp, const char *attrib, bool *p, size_t size)
{
	GA_RWHandleI handle = gdp.findPointAttribute(attrib);
	if (gdp.getNumPoints() != size)
		RESIZE(gdp, size);
	for (size_t i = 0; i < size; i++)
		handle.set(gdp.pointOffset(i), p[i]);
}
inline void SYNC_I(GU_Detail &gdp, const char *attrib, int *p, size_t size)
{
	GA_RWHandleI handle = gdp.findPointAttribute(attrib);
	if (gdp.getNumPoints() != size)
		RESIZE(gdp, size);
	for (size_t i = 0; i < size; i++)
		handle.set(gdp.pointOffset(i), p[i]);
}
inline void SYNC_F(GU_Detail &gdp, const char *attrib, float *p, size_t size)
{
	GA_RWHandleF handle = gdp.findPointAttribute(attrib);
	if (gdp.getNumPoints() != size)
		RESIZE(gdp, size);
	for (size_t i = 0; i < size; i++)
		handle.set(gdp.pointOffset(i), p[i]);
}
inline void SYNC_V3(GU_Detail &gdp, const char *attrib, float *p, size_t size)
{
	GA_RWHandleV3 handle = gdp.findPointAttribute(attrib);
	if (gdp.getNumPoints() != size)
		RESIZE(gdp, size);
	for (size_t i = 0; i < size; i++)
		handle.set(gdp.pointOffset(i), UT_Vector3(p[3 * i], p[3 * i + 1], p[3 * i + 2]));
}

#endif //DFSPH_UTILS_H
