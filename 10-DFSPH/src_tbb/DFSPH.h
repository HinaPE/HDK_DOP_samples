#ifndef DFSPH_DFSPH_TBB_H
#define DFSPH_DFSPH_TBB_H

#include <vector>
#include <memory>
#include <cmath>
#include <UT/UT_Vector3.h>

class GU_Detail;

namespace HinaPE::TBB
{
using float3 = UT_Vector3;
using ScalarArray = std::vector<float>;
using IntegerArray = std::vector<float>;
using Vector3Array = std::vector<float3>;
struct FluidData
{
	Vector3Array x;
	Vector3Array v;
	Vector3Array a;
	ScalarArray m;
	ScalarArray V;
	ScalarArray rho;

	// temp buffers
	ScalarArray factor;
	ScalarArray density_adv;
	ScalarArray nn;
	ScalarArray tmp;
};

struct DFSPHParams
{
	size_t MAX_ITERATIONS = 100;
	float REST_DENSITY = 1000.f;
	float PARTICLE_RADIUS = 0.01f;
	float SURFACE_TENSION = 0.01f;
	float VISCOSITY = 0.01f;
	float DEFAULT_V = std::powf(2 * PARTICLE_RADIUS, 3);
	float DEFAULT_M = REST_DENSITY * DEFAULT_V;

	float3 GRAVITY = {0, -9.8f, 0};
	float3 MaxBound;
	bool TOP_OPEN = true;
};

struct DFSPH : DFSPHParams
{
	DFSPH(float _kernel_radius);
	void solve(float dt, const GU_Detail *gdp);
	std::shared_ptr<FluidData> Fluid;
private:
	float size;
	float kr;
};
}

#endif //DFSPH_DFSPH_TBB_H
