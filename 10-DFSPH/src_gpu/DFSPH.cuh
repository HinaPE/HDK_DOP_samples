#ifndef DFSPH_DFSPH_CUDA_CUH
#define DFSPH_DFSPH_CUDA_CUH

#include <vector>
#include <memory>
#include <cmath>
#include <thrust/universal_vector.h>

namespace cuNSearch { class NeighborhoodSearch; }

namespace HinaPE::CUDA
{
using ScalarArray = thrust::universal_vector<float>;
using Vector3Array = thrust::universal_vector<float3>;
struct FluidGPU
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
	// set outside
	float KERNEL_RADIUS = 0.04;
	float SURFACE_TENSION = 0.01f;
	float VISCOSITY = 0.01f;
	float3 MaxBound = {1, 1, 1};
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
	void set_gpu_constants();
	void resize(size_t n);
	void solve(float dt);

	std::shared_ptr<FluidGPU> Fluid;
	size_t size;

private:
	std::shared_ptr<cuNSearch::NeighborhoodSearch> Searcher;
	size_t fluid_idx;
	bool need_reload;
};
}

#endif //DFSPH_DFSPH_CUDA_CUH
