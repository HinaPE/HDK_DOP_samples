#include "PBF.h"
#include <GU/GU_Detail.h>
#include <GU/GU_NeighbourList.h>
#include <UT/UT_ParallelUtil.h>
#include <numeric>

// @formatter:off
namespace HinaPE::TBB
{
struct Cubic
{
	inline static void SetRadius(float r) { _r = r; const float pi = 3.14159265358979323846f; const float h3 = _r * _r * _r; _k = 8.f / (pi * h3); _l = 48.f / (pi * h3); _W_0 = W(float3{0, 0, 0}); }
	inline static float W(const float r) { float res = 0.f; const float q = r / _r; if (q <= 1.0) { if (q <= 0.5) { const float q2 = q * q; const float q3 = q2 * q; res = _k * (6.f * q3 - 6.f * q2 + 1.f); } else { res = _k * (2.f * pow(1.f - q, 3.f)); } } return res; }
	inline static float W(const float3 &r) { return W(r.length()); }
	inline static float3 gradW(const float3 &r) { float3 res; const float rl = r.length(); const float q = rl / _r; if ((rl > 1.0e-9) && (q <= 1.f)) { float3 gradq = r / rl; gradq /= _r; if (q <= 0.5) { res = _l * q * (3.f * q - 2.f) * gradq; } else { const float factor = 1.f - q; res = _l * (-factor * factor) * gradq; } } else res = float3{0, 0, 0}; return res; }
	inline static float W_zero() { return _W_0; }
	inline static float _r;
	inline static float _k;
	inline static float _l;
	inline static float _W_0;
};
using Kernel = Cubic;
}
// @formatter:on

inline void parallel_for(size_t n, const std::function<void(size_t)> &f)
{
	UTparallelForEachNumber((int) n, [&](const UT_BlockedRange<int> &range)
	{
		for (size_t i = range.begin(); i != range.end(); ++i) { f(i); }
	});
}

HinaPE::TBB::PBF::PBF() : size(0)
{
	Fluid = std::make_shared<FluidData>();
	Kernel::SetRadius(0.04);
}
void HinaPE::TBB::PBF::solve(float dt, const GU_Detail *gdp)
{
	if (Kernel::_r != KERNEL_RADIUS)
		Kernel::SetRadius(KERNEL_RADIUS);
}

#ifdef TEST_PBF
#include <iostream>
int main()
{
	std::cout << "Hello, TBB!" << std::endl;
	return 0;
}
#endif
