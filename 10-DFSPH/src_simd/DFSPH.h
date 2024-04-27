#ifndef DFSPH_DFSPH_SIMD_H
#define DFSPH_DFSPH_SIMD_H

#include <vector>
#include <array>
#include <memory>

namespace HinaPE::SIMD
{
using ScalarArray = std::vector<float>;
using Vector3Array = std::vector<std::array<float, 3>>;

struct FluidSIMD
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

	std::shared_ptr<FluidSIMD> Fluid;
private:
	size_t size;
	float kernel_radius;
	size_t fluid_idx;
};
}

#endif //DFSPH_DFSPH_SIMD_H
