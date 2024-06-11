#include "GAS_Smoke_Solver.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_ForceGravity.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_FieldUtils.h>
#include <GAS/GAS_ProjectNonDivergent.h>
#include <GAS/GAS_Diffuse.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>

#include "src/advect.h"
#include "src/diffuse.h"
#include "src/pressure.h"
#include "src/source.h"

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_SOURCE static PRM_Name SourceName(GAS_NAME_SOURCE, "Source"); static PRM_Default SourceNameDefault(0, GAS_NAME_SOURCE); PRMs.emplace_back(PRM_STRING, 1, &SourceName, &SourceNameDefault);
#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);
#define ACTIVATE_GAS_TEMPERATURE static PRM_Name TemperatureName(GAS_NAME_TEMPERATURE, "Temperature"); static PRM_Default TemperatureNameDefault(0, GAS_NAME_TEMPERATURE); PRMs.emplace_back(PRM_STRING, 1, &TemperatureName, &TemperatureNameDefault);
#define ACTIVATE_GAS_COLLISION static PRM_Name CollisionName(GAS_NAME_COLLISION, "Collision"); static PRM_Default CollisionNameDefault(0, GAS_NAME_COLLISION); PRMs.emplace_back(PRM_STRING, 1, &CollisionName, &CollisionNameDefault);
#define ACTIVATE_GAS_VELOCITY static PRM_Name VelocityName(GAS_NAME_VELOCITY, "Velocity"); static PRM_Default VelocityNameDefault(0, GAS_NAME_VELOCITY); PRMs.emplace_back(PRM_STRING, 1, &VelocityName, &VelocityNameDefault);
#define ACTIVATE_GAS_PRESSURE static PRM_Name PressureName(GAS_NAME_PRESSURE, "Pressure"); static PRM_Default PressureNameDefault(0, GAS_NAME_PRESSURE); PRMs.emplace_back(PRM_STRING, 1, &PressureName, &PressureNameDefault);
#define ACTIVATE_GAS_DIVERGENCE static PRM_Name DivergenceName(GAS_NAME_DIVERGENCE, "Divergence"); static PRM_Default DivergenceNameDefault(0, GAS_NAME_DIVERGENCE); PRMs.emplace_back(PRM_STRING, 1, &DivergenceName, &DivergenceNameDefault);

#define GAS_NAME_VELOCITY_SWAP        "velocity_swap"
#define ACTIVATE_GAS_VELOCITY_SWAP static PRM_Name VelocitySwapName(GAS_NAME_VELOCITY_SWAP, "VelocitySwap"); static PRM_Default VelocitySwapNameDefault(0, GAS_NAME_VELOCITY_SWAP); PRMs.emplace_back(PRM_STRING, 1, &VelocitySwapName, &VelocitySwapNameDefault);

#define GAS_NAME_PRESSURE_SWAP        "pressure_swap"
#define ACTIVATE_GAS_PRESSURE_SWAP static PRM_Name PressureSwapName(GAS_NAME_PRESSURE_SWAP, "PressureSwap"); static PRM_Default PressureSwapNameDefault(0, GAS_NAME_PRESSURE_SWAP); PRMs.emplace_back(PRM_STRING, 1, &PressureSwapName, &PressureSwapNameDefault);

#define GAS_NAME_DIVERGENCE_SWAP        "divergence_swap"
#define ACTIVATE_GAS_DIVERGENCE_SWAP static PRM_Name DivergenceSwapName(GAS_NAME_DIVERGENCE_SWAP, "DivergenceSwap"); static PRM_Default DivergenceSwapNameDefault(0, GAS_NAME_DIVERGENCE_SWAP); PRMs.emplace_back(PRM_STRING, 1, &DivergenceSwapName, &DivergenceSwapNameDefault);

#define PARAMETER_FLOAT(NAME, DEFAULT_VALUE) static PRM_Name NAME(#NAME, #NAME);static PRM_Default Default##NAME(DEFAULT_VALUE);PRMs.emplace_back(PRM_FLT, 1, &NAME, &Default##NAME);

void GAS_Smoke_Solver::initializeSubclass()
{
	SIM_Data::initializeSubclass();
	this->Advection = std::make_shared<HinaPE::AdvectionSolver>();
	this->Diffusion = std::make_shared<HinaPE::DiffusionSolver>();
	this->Poisson = std::make_shared<HinaPE::PoissonSolver>();
	this->Source = std::make_shared<HinaPE::SourceSolver>();
	this->D_Swap = std::make_shared<SIM_RawField>();
	this->T_Swap = std::make_shared<SIM_RawField>();
}
void GAS_Smoke_Solver::makeEqualSubclass(const SIM_Data *source)
{
	SIM_Data::makeEqualSubclass(source);
	this->Advection = ((GAS_Smoke_Solver *) source)->Advection;
	this->Diffusion = ((GAS_Smoke_Solver *) source)->Diffusion;
	this->Poisson = ((GAS_Smoke_Solver *) source)->Poisson;
	this->Source = ((GAS_Smoke_Solver *) source)->Source;
	this->D_Swap = ((GAS_Smoke_Solver *) source)->D_Swap;
	this->T_Swap = ((GAS_Smoke_Solver *) source)->T_Swap;
}
const SIM_DopDescription *GAS_Smoke_Solver::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	ACTIVATE_GAS_SOURCE
	ACTIVATE_GAS_DENSITY
	ACTIVATE_GAS_TEMPERATURE
	ACTIVATE_GAS_COLLISION
	ACTIVATE_GAS_VELOCITY
	ACTIVATE_GAS_PRESSURE
	ACTIVATE_GAS_DIVERGENCE
	ACTIVATE_GAS_VELOCITY_SWAP
	ACTIVATE_GAS_PRESSURE_SWAP
	ACTIVATE_GAS_DIVERGENCE_SWAP
	PARAMETER_FLOAT(Viscosity, 0.01)
	PARAMETER_FLOAT(Diffusion, 0.01)
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
bool GAS_Smoke_Solver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_SOURCE);
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);
	SIM_ScalarField *T = getScalarField(obj, GAS_NAME_TEMPERATURE);
	SIM_VectorField *V = getVectorField(obj, GAS_NAME_VELOCITY);
	SIM_ScalarField *P = getScalarField(obj, GAS_NAME_PRESSURE);
	SIM_ScalarField *Div = getScalarField(obj, GAS_NAME_DIVERGENCE);

	SIM_VectorField *V_Swap = getVectorField(obj, GAS_NAME_VELOCITY_SWAP);
	SIM_ScalarField *P_Swap = getScalarField(obj, GAS_NAME_PRESSURE_SWAP);
	SIM_ScalarField *Div_Swap = getScalarField(obj, GAS_NAME_DIVERGENCE_SWAP);

	if (!G || !S || !D || !T || !V || !P || !Div || !V_Swap || !P_Swap || !Div_Swap)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	float visc = getViscosity();
	float diff = getDiffusion();
	float h = D->getField()->getVoxelSize().x();
	float factor = timestep * diff / (h * h);

//	V->buildDivergenceCenter(*Div->getField());
//	for (int axis: {0, 1, 2})
//	{
//		V_Swap->getField(axis)->fieldNC()->copyData(*V->getField(axis)->field());
//	}
//	V_Swap->projectToNonDivergent();
//	V_Swap->buildDivergenceCenter(*Div_Swap->getField());

	Source->Buoyancy(V, D->getField(), T->getField(), T->getField()->average(), timestep);
	V->projectToNonDivergent();
	V->enforceBoundary();
	V->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
	V->enforceBoundary();
	V->projectToNonDivergent();
	V->enforceBoundary();

	Source->Emit(D->getField()->fieldNC(), S->getField()->field(), V, UT_Vector3(0, 0, 0));
	Source->Emit(T->getField()->fieldNC(), S->getField()->field(), V, UT_Vector3(0, 0, 0));
	SIM_RawIndexField MARKER;
	MARKER.init(D->getVoxelSample(), D->getOrig(), D->getSize(), D->getField()->field()->getXRes(), D->getField()->field()->getYRes(), D->getField()->field()->getZRes());
	MARKER.fieldNC()->constant(1);
	SIM_RawField D_Swap1(*D->getField());
	Diffusion->PCG(D->getField()->fieldNC(), D_Swap1.field(), MARKER.field(), factor);
	SIM_RawField T_Swap1(*T->getField());
	Diffusion->PCG(T->getField()->fieldNC(), T_Swap1.field(), MARKER.field(), factor);
	D->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
	D->enforceBoundary();
	T->advect(V, -timestep, nullptr, SIM_ADVECT_MIDPOINT, 1.0f);
	T->enforceBoundary();

	return true;
}
