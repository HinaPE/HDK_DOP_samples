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
};

struct DFSPHParams
{
	// set outside
	float KERNEL_RADIUS = 0.04;
	float SURFACE_TENSION = 0.01f;
	float VISCOSITY = 0.01f;
	float3 MaxBound;
	bool TOP_OPEN = true;
	size_t DIVERGENCE_ITERS = 0;
	size_t PRESSURE_ITERS = 0;

	// immutable
	const float PARTICLE_RADIUS = 0.01f;
	const float3 GRAVITY = {0, -9.8f, 0};
	const size_t MAX_ITERATIONS = 100;
	const float REST_DENSITY = 1000.f;
	const float DEFAULT_V = std::powf(2 * PARTICLE_RADIUS, 3);
	const float DEFAULT_M = REST_DENSITY * DEFAULT_V;
};

struct DFSPH : DFSPHParams
{
	DFSPH();
	void solve(float dt, const GU_Detail *gdp);
	std::shared_ptr<FluidData> Fluid;

private:
	float size;
};
}

#endif //DFSPH_DFSPH_TBB_H
