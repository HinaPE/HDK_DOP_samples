#ifndef DFSPH_DFSPH_TBB_H
#define DFSPH_DFSPH_TBB_H

#include <vector>
#include <memory>

namespace HinaPE::TBB
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

};
}

#endif //DFSPH_DFSPH_TBB_H
