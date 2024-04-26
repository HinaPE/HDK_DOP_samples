#include "GAS_DFSPH_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>

#include "SIM_DFSPH_Particles.h"
#include "src_tbb/DFSPH.h"
#include "src_simd/DFSPH.h"
#include "src_gpu/DFSPH.cuh"

void GAS_DFSPH_Solver::initializeSubclass() { SIM_Data::initializeSubclass(); }
void GAS_DFSPH_Solver::makeEqualSubclass(const SIM_Data *source) { SIM_Data::makeEqualSubclass(source); }
const SIM_DopDescription *GAS_DFSPH_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	TARGET_SOLVE_GEOMETRY(SIM_DFSPH_Particles)
	PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   ENGLISH_NAME,
								   COMMON_NAME,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	setGasDescription(DESC);
	return &DESC;
}
bool GAS_DFSPH_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_DFSPH_Particles *particles = SIM_DATA_GET(*obj, SIM_DFSPH_Particles::DATANAME, SIM_DFSPH_Particles);
	if (!particles)
		return false;

	SIM_GeometryAutoWriteLock lock(particles);
	GU_Detail &gdp = lock.getGdp();
//	SYNC_V3(gdp, SIM_ATTRIBUTE_NAME_P, nullptr, 2);

	GA_Offset pt_off;
	GA_FOR_ALL_PTOFF(&gdp, pt_off)
		{
			auto pos = gdp.getPos3(pt_off);
			gdp.setPos3(pt_off, pos + UT_Vector3(0, -0.1, 0));
		}

	return true;
}
