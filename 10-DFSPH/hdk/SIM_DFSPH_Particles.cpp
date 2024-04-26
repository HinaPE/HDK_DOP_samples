#include "SIM_DFSPH_Particles.h"

#include <SIM/SIM_DopDescription.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

void SIM_DFSPH_Particles::initializeSubclass() { SIM_GeometryCopy::initializeSubclass(); }
void SIM_DFSPH_Particles::makeEqualSubclass(const SIM_Data *source) { SIM_GeometryCopy::makeEqualSubclass(source); }

const SIM_DopDescription *SIM_DFSPH_Particles::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	PARAMETER_BOOL(BBB, 1)
	PARAMETER_INT(III, 1)
	PARAMETER_FLOAT(FFF, 1.0)
	PARAMETER_STRING(STR, "string")
	PARAMETER_VECTOR_INT_N(VEC3I, 3, 1, 2, 3)
	PARAMETER_VECTOR_FLOAT_N(VEC3F, 3, 1.0, 2.0, 3.0)
//	PARAMETER_POSITION
//	PARAMETER_ROTATION
	PARAMETER_SUB_POSITION_PATH
	PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   ENGLISH_NAME,
								   COMMON_NAME,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	return &DESC;
}

GU_ConstDetailHandle SIM_DFSPH_Particles::getGeometrySubclass() const
{
	if (_handle.isNull())
	{
		GU_Detail *gdp = new GU_Detail();
		GEOMETRY_POINT_ATTRIBUTE_INT("INT");
		GEOMETRY_POINT_ATTRIBUTE_FLOAT("FLOAT");
		GEOMETRY_POINT_ATTRIBUTE_STRING("STRING");
		GEOMETRY_POINT_ATTRIBUTE_VECTOR_FLOAT_3("VECTOR3");
		_handle.allocateAndSet(gdp);
	}
	return _handle;
}
