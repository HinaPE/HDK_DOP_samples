#include <UT/UT_DSOVersion.h> // Very Important!!! Include this first

#include "SIM_DFSPH_Particles.h"
#include "GAS_DFSPH_Solver.h"

void initializeSIM(void *)
{
	IMPLEMENT_DATAFACTORY(SIM_DFSPH_Particles)
	IMPLEMENT_DATAFACTORY(GAS_DFSPH_Solver)
}
