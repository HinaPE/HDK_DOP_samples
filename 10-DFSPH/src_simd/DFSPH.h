#ifndef DFSPH_DFSPH_SIMD_H
#define DFSPH_DFSPH_SIMD_H

#include <vector>
#include <memory>

namespace HinaPE::SIMD
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
};

using float3 = std::array<float, 3>;
using ScalarArrayGPU = std::vector<float>;
using Vector3ArrayGPU = std::vector<float3>;
using FluidGPU = IFluid<ScalarArrayGPU, Vector3ArrayGPU>;

struct DFSPH
{
	DFSPH(float _kernel_radius);
	void resize(size_t n);
	void solve(float dt);

private:
	std::shared_ptr<FluidGPU> Fluid;
	size_t size;
	float kernel_radius;
	size_t fluid_idx;
};
}

#endif //DFSPH_DFSPH_SIMD_H
