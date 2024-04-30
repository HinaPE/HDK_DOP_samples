#include "GAS_Rigid_Apply_Impulse.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#include <RBD/RBD_Object.h>
#include <RBD/RBD_Solver.h>

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);
#define PARAMETER_VECTOR_FLOAT_N(NAME, SIZE, ...) static PRM_Name NAME(#NAME, #NAME); static std::array<PRM_Default, SIZE> Default##NAME{__VA_ARGS__}; PRMs.emplace_back(PRM_FLT, SIZE, &NAME, Default##NAME.data());

const SIM_DopDescription *GAS_Rigid_Apply_Impulse::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	PARAMETER_FLOAT(Force, 1.f)
	PARAMETER_VECTOR_FLOAT_N(Position, 3, 0.f, 0.f, 0.f)
	PARAMETER_VECTOR_FLOAT_N(Direction, 3, 0.f, 1.f, 0.f)
	PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   DOP_NAME,
								   DOP_ENGLISH,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	setGasDescription(DESC);
	return &DESC;
}
bool GAS_Rigid_Apply_Impulse::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_ObjectArray affectors;
	obj->getAffectors(affectors, "SIM_RelationshipCollide");

	exint num_affectors = affectors.entries();
	for (int i = 0; i < num_affectors; ++i)
	{
		SIM_Object *rbd_obj = affectors(i);
		if (rbd_obj->getName().equal(obj->getName()))
			continue;
		const RBD_Solver *rbd_solver = SIM_DATA_CASTCONST(rbd_obj->getSolver(), RBD_Solver);
		if (rbd_solver)
		{
			RBD_Object RBD(rbd_solver, rbd_obj);
			RBD.applyImpulse(getPosition(), getForce(), getDirection());
		}
	}
	return true;
}
