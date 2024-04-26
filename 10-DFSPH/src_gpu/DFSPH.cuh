#ifndef DFSPH_DFSPH_CUDA_CUH
#define DFSPH_DFSPH_CUDA_CUH

#include <thrust/universal_vector.h>

namespace cuNSearch { class NeighborhoodSearch; }

namespace HinaPE::CUDA
{
template<typename ScalarArray, typename Vector3Array>
struct IFluid
{
	Vector3Array x;
	Vector3Array v;
	Vector3Array a;
	ScalarArray m;
	ScalarArray V;
	ScalarArray rho;

	// temp
	ScalarArray factor;
	ScalarArray density_adv;
};

using ScalarArrayGPU = thrust::universal_vector<float>;
using Vector3ArrayGPU = thrust::universal_vector<float3>;
using FluidGPU = IFluid<ScalarArrayGPU, Vector3ArrayGPU>;

struct DFSPH
{
	DFSPH(float _kernel_radius);
	void resize(size_t n);
	void solve(float dt);

private:
	std::shared_ptr<FluidGPU> Fluid;
	std::shared_ptr<cuNSearch::NeighborhoodSearch> Searcher;
	size_t size;
	size_t fluid_idx;
};
}

#endif //DFSPH_DFSPH_CUDA_CUH
