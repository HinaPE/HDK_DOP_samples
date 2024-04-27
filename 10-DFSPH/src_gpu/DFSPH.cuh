#ifndef DFSPH_DFSPH_CUDA_CUH
#define DFSPH_DFSPH_CUDA_CUH

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

struct DFSPH
{
	DFSPH(float _kernel_radius);
	void resize(size_t n);
	void solve(float dt);
	void solve_test(float dt);

	std::shared_ptr<FluidGPU> Fluid;

private:
	std::shared_ptr<cuNSearch::NeighborhoodSearch> Searcher;
	size_t size;
	size_t fluid_idx;
	bool need_reload;
};
}

#endif //DFSPH_DFSPH_CUDA_CUH
