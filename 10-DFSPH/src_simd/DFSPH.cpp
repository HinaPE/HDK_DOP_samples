#include "DFSPH.h"

#include <GU/GU_Detail.h>
#include <GU/GU_NeighbourList.h>
#include <UT/UT_ParallelUtil.h>
#include <VM/VM_Math.h>

// @formatter:off
namespace HinaPE::SIMD
{
struct Cubic
{
	inline static void SetRadius(float r) { _r = r; const float pi = 3.14159265358979323846f; const float h3 = _r * _r * _r; _k = 8.f / (pi * h3); _l = 48.f / (pi * h3); _W_0 = W(float3{0, 0, 0}); }
	inline static float W(const float r) { float res = 0.0; const float q = r / _r; if (q <= 1.0) { if (q <= 0.5) { const float q2 = q * q; const float q3 = q2 * q; res = _k * (6.f * q3 - 6.f * q2 + 1.f); } else { res = _k * (2.f * pow(1.f - q, 3.f)); } } return res; }
	inline static float W(const float3 &r) { return W(r.length()); }
	inline static float3 gradW(const float3 &r) { float3 res; const float rl = r.length(); const float q = rl / _r; if ((rl > 1.0e-9) && (q <= 1.0)) { float3 gradq = r / rl; gradq /= _r; if (q <= 0.5) { res = _l * q * (3.f * q - 2.f) * gradq; } else { const float factor = 1.f - q; res = _l * (-factor * factor) * gradq; } } else res = float3{0, 0, 0}; return res; }
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

HinaPE::SIMD::DFSPH::DFSPH()
{
	Fluid = std::make_shared<FluidData>();
	Kernel::SetRadius(0.04);
}
void HinaPE::SIMD::DFSPH::solve(float dt, const GU_Detail *gdp)
{
	VM_Math::forceSIMD(true);

	if (Kernel::_r != KERNEL_RADIUS)
		Kernel::SetRadius(KERNEL_RADIUS);

	// ==================== 1. Refresh Fluid Data ====================
	static std::vector<UT_Array<GA_Offset>> NL;
	if (size < gdp->getNumPoints())
	{
		auto new_size = gdp->getNumPoints();
		Fluid->x.resize(3 * new_size, 0);
		Fluid->v.resize(3 * new_size, 0);
		Fluid->a.resize(3 * new_size, 0);
		Fluid->m.resize(new_size, DEFAULT_M);
		Fluid->V.resize(new_size, DEFAULT_V);
		Fluid->rho.resize(new_size, 0);
		Fluid->factor.resize(new_size, 0);
		Fluid->density_adv.resize(new_size, 0);
		Fluid->nn.resize(new_size, 0);
		Fluid->tmpi.resize(new_size, 0);
		Fluid->tmpv.resize(3 * new_size, 0);
		Fluid->tmp.resize(new_size, 0);
		NL.resize(new_size);
		Fluid->nn.resize(new_size, 0);
		size = new_size;

		{
			GA_Offset pt_off;
			GA_FOR_ALL_PTOFF(gdp, pt_off)
				{
					GA_Index pt_idx = gdp->pointIndex(pt_off);
					UT_Vector3 pos = gdp->getPos3(pt_off);
					Fluid->x[3 * pt_idx + 0] = pos.x();
					Fluid->x[3 * pt_idx + 1] = pos.y();
					Fluid->x[3 * pt_idx + 2] = pos.z();
				}
		}
	}
	size_t remains = size % PATCHES;



	// ==================== 2. Build Neighbors ====================
	GU_NeighbourList Searcher;
	GU_NeighbourListParms NLP;
	NLP.setRadius(KERNEL_RADIUS);
	NLP.setOverrideRadius(true);
	NLP.setMode(GU_NeighbourListParms::InteractionMode::UNIFORM);
	Searcher.build(gdp, NLP);
	parallel_for(size, [&](size_t i)
	{
		UT_Array<GA_Offset> neighbor_list;
		Searcher.getNeighbours(i, gdp, neighbor_list);
		NL[i] = neighbor_list;
		Fluid->nn[i] = neighbor_list.size();
	});



	// ==================== 3. Compute Density, Factor ====================
	std::fill(Fluid->tmp.begin(), Fluid->tmp.end(), 0); // sum_grad
	std::fill(Fluid->tmpv.begin(), Fluid->tmpv.end(), 0); // grad_i
	for (int i = 0; i < size - remains; i += PATCHES)
	{
		VM_Math::mul(&Fluid->rho[i + 0], &Fluid->V[i + 0], Kernel::W_zero(), PATCHES);

		for (int patch = 0; patch < PATCHES; ++patch)
		{
			float3 x_i{Fluid->x[3 * (i + patch) + 0], Fluid->x[3 * (i + patch) + 1], Fluid->x[3 * (i + patch) + 2]};
			for (auto j: NL[i + patch])
			{
				float3 x_j{Fluid->x[3 * j + 0], Fluid->x[3 * j + 1], Fluid->x[3 * j + 2]};
				Fluid->rho[i + patch] += Fluid->V[j] * Kernel::W(x_i - x_j);

				float3 grad_j = -Fluid->V[j] * Kernel::gradW(x_i - x_j);
				Fluid->tmp[i + patch] += grad_j.length2();
				Fluid->tmpv[3 * (i + patch) + 0] -= grad_j.x();
				Fluid->tmpv[3 * (i + patch) + 1] -= grad_j.y();
				Fluid->tmpv[3 * (i + patch) + 2] -= grad_j.z();
			}
			float3 grad_i_patch = {Fluid->tmpv[3 * (i + patch) + 0], Fluid->tmpv[3 * (i + patch) + 1], Fluid->tmpv[3 * (i + patch) + 2]};
			Fluid->tmp[i + patch] += grad_i_patch.length2();
		}
		VM_Math::mulSIMD(&Fluid->rho[i + 0], &Fluid->rho[i + 0], REST_DENSITY, PATCHES);
		VM_Math::setSIMD(&Fluid->factor[i + 0], -1.f, PATCHES);
		VM_Math::safedivSIMD(&Fluid->factor[i + 0], &Fluid->factor[i + 0], &Fluid->tmp[i + 0], PATCHES);
	}
	for (int i = size - remains; i < size; ++i)
	{
		float3 x_i{Fluid->x[3 * i + 0], Fluid->x[3 * i + 1], Fluid->x[3 * i + 2]};

		float sum_grad = 0.;
		float3 grad_i{0, 0, 0};
		Fluid->rho[i] = Fluid->V[i] * Kernel::W_zero();
		for (auto j: NL[i])
		{
			float3 x_j{Fluid->x[3 * j + 0], Fluid->x[3 * j + 1], Fluid->x[3 * j + 2]};
			Fluid->rho[i] += Fluid->V[j] * Kernel::W(x_i - x_j);
			float3 grad_j = -Fluid->V[j] * Kernel::gradW(x_i - x_j);
			sum_grad += grad_j.length2();
			grad_i -= grad_j;
		}
		Fluid->rho[i] *= REST_DENSITY;
		sum_grad += grad_i.length2();
		if (sum_grad > 1e-6)
			Fluid->factor[i] = -1.f / sum_grad;
		else
			Fluid->factor[i] = 0;
	}
	std::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(), [](float f) { return f == -1 ? 0 : f; });



	// ==================== 5. Non-Pressure Force and Predict Velocity ====================
	constexpr float d = 10;
	const float diameter = 2 * PARTICLE_RADIUS;
	const float diameter2 = diameter * diameter;
	const float kr2 = KERNEL_RADIUS * KERNEL_RADIUS;
	parallel_for(size, [&](size_t i)
	{
		float3 dv = GRAVITY;
		for (auto j: NL[i])
		{
			// Surface Tension
			const float3 r{Fluid->x[3 * i + 0] - Fluid->x[3 * j + 0], Fluid->x[3 * i + 1] - Fluid->x[3 * j + 1], Fluid->x[3 * i + 2] - Fluid->x[3 * j + 2]};
			const float r2 = r.dot(r);
			const float r1 = std::sqrt(r2);
			if (r2 > diameter2)
				dv -= SURFACE_TENSION / Fluid->m[i] * Fluid->m[j] * r * Kernel::W(r1);
			else
				dv -= SURFACE_TENSION / Fluid->m[i] * Fluid->m[j] * r * Kernel::W(diameter);

			// Fluid Viscosity
			float v_xy = (Fluid->v[3 * i + 0] - Fluid->v[3 * j + 0]) * r.x() + (Fluid->v[3 * i + 1] - Fluid->v[3 * j + 1]) * r.y() + (Fluid->v[3 * i + 2] - Fluid->v[3 * j + 2]) * r.z();
			float3 f_v = d * VISCOSITY * (Fluid->m[j] / (Fluid->rho[j])) * v_xy / (r2 + 0.01f * kr2) * Kernel::gradW(r);
			dv += f_v;
		}
		Fluid->a[3 * i + 0] = dv.x();
		Fluid->a[3 * i + 1] = dv.y();
		Fluid->a[3 * i + 2] = dv.z();
	});
	for (int i = 0; i < size - remains; i += PATCHES)
	{
		VM_Math::mulSIMD(&Fluid->tmpv[3 * i + 0], &Fluid->a[3 * i + 0], dt, PATCHES * 3);
		VM_Math::addSIMD(&Fluid->v[3 * i + 0], &Fluid->v[3 * i + 0], &Fluid->tmpv[3 * i + 0], PATCHES * 3);
	}
	for (int i = size - remains; i < size; ++i)
	{
		Fluid->v[3 * i + 0] += dt * Fluid->a[3 * i + 0];
		Fluid->v[3 * i + 1] += dt * Fluid->a[3 * i + 1];
		Fluid->v[3 * i + 2] += dt * Fluid->a[3 * i + 2];
	}


	// ==================== 7. Advect ====================
	for (int i = 0; i < size - remains; i += PATCHES)
	{
		VM_Math::mulSIMD(&Fluid->tmpv[3 * i + 0], &Fluid->v[3 * i + 0], dt, PATCHES * 3);
		VM_Math::addSIMD(&Fluid->x[3 * i + 0], &Fluid->x[3 * i + 0], &Fluid->tmpv[3 * i + 0], PATCHES * 3);
	}
	for (int i = size - remains; i < size; ++i)
	{
		Fluid->x[3 * i + 0] += dt * Fluid->v[3 * i + 0];
		Fluid->x[3 * i + 1] += dt * Fluid->v[3 * i + 1];
		Fluid->x[3 * i + 2] += dt * Fluid->v[3 * i + 2];
	}


	// ==================== 8. Enforce Boundary ====================
	parallel_for(size, [&](size_t i)
	{
		float3 normal{0, 0, 0};
		if (Fluid->x[3 * i + 0] > MaxBound.x())
		{
			Fluid->x[3 * i + 0] = MaxBound.x();
			normal.x() += 1;
		}
		if (Fluid->x[3 * i + 0] < -MaxBound.x())
		{
			Fluid->x[3 * i + 0] = -MaxBound.x();
			normal.x() -= 1;
		}
		if (!TOP_OPEN)
		{
			if (Fluid->x[3 * i + 1] > MaxBound.y())
			{
				Fluid->x[3 * i + 1] = MaxBound.y();
				normal.y() += 1;
			}
		}
		if (Fluid->x[3 * i + 1] < -MaxBound.y())
		{
			Fluid->x[3 * i + 1] = -MaxBound.y();
			normal.y() -= 1;
		}
		if (Fluid->x[3 * i + 2] > MaxBound.z())
		{
			Fluid->x[3 * i + 2] = MaxBound.z();
			normal.z() += 1;
		}
		if (Fluid->x[3 * i + 2] < -MaxBound.z())
		{
			Fluid->x[3 * i + 2] = -MaxBound.z();
			normal.z() -= 1;
		}
		normal.normalize();
		constexpr float c_f = 0.5f;
		float3 v_i{Fluid->v[3 * i + 0], Fluid->v[3 * i + 1], Fluid->v[3 * i + 2]};
		float3 res = (1.f + c_f) * v_i.dot(normal) * normal;
		Fluid->v[3 * i + 0] -= res.x();
		Fluid->v[3 * i + 1] -= res.y();
		Fluid->v[3 * i + 2] -= res.z();
	});
}

//void HinaPE::SIMD::DFSPH::solve(float dt, GU_Detail &gdp)
//{
//	// ==================== 1. Build Neighbors ====================
//	Searcher->build(&gdp, *nlp);
//
//
//
//	// ==================== 2. Compute Density and Factor ====================
//	parallel_for(size, [&](size_t i)
//	{
//		UT_Array<GA_Offset> neighbor_list;
//		Searcher->getNeighbours(i, &gdp, neighbor_list);
//		Fluid->nn[i] = neighbor_list.size();
//	});
//
//	std::transform(Fluid->a.begin(), Fluid->a.end(), Fluid->a.begin(), [](std::array<float, 3> a) { return std::array<float, 3>{0, -9.8, 0}; });
//
////	std::chrono::duration<double, std::milli> duration1, duration2, duration3, duration4;
////	{
////		auto start = std::chrono::high_resolution_clock::now();
////		size_t patch = size / 16;
////		size_t left = size % 16;
////
////		parallel_for(patch, [&](size_t p)
////		{
////			size_t i = p * 16;
////			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->a[i][0], dt, 16 * 3);
////			VM_Math::addSIMD(&Fluid->v[i][0], &Fluid->v[i][0], &Fluid->tmpv[i][0], 16 * 3);
////			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->v[i][0], dt, 16 * 3);
////			VM_Math::addSIMD(&Fluid->x[i][0], &Fluid->x[i][0], &Fluid->tmpv[i][0], 16 * 3);
////		});
////
////		for (size_t i = size - left; i < size; ++i)
////		{
////			Fluid->v[i][0] += dt * Fluid->a[i][0];
////			Fluid->v[i][1] += dt * Fluid->a[i][1];
////			Fluid->v[i][2] += dt * Fluid->a[i][2];
////
////			Fluid->x[i][0] += dt * Fluid->v[i][0];
////			Fluid->x[i][1] += dt * Fluid->v[i][1];
////			Fluid->x[i][2] += dt * Fluid->v[i][2];
////		}
////
////		auto end = std::chrono::high_resolution_clock::now();
////		duration1 = end - start;
////		std::cout << "parallel + SIMD Time: " << duration1.count() << " ms" << std::endl;
////	}
////
////	{
////		auto start = std::chrono::high_resolution_clock::now();
////		parallel_for(size, [&](size_t i)
////		{
////			Fluid->v[i][0] += dt * Fluid->a[i][0];
////			Fluid->v[i][1] += dt * Fluid->a[i][1];
////			Fluid->v[i][2] += dt * Fluid->a[i][2];
////
////			Fluid->x[i][0] += dt * Fluid->v[i][0];
////			Fluid->x[i][1] += dt * Fluid->v[i][1];
////			Fluid->x[i][2] += dt * Fluid->v[i][2];
////		});
////		auto end = std::chrono::high_resolution_clock::now();
////		duration2 = end - start;
////		std::cout << "parallel Time: " << duration2.count() << " ms" << std::endl;
////	}
////
////	{
////		auto start = std::chrono::high_resolution_clock::now();
////		for (int i = 0; i < size; ++i)
////		{
////			Fluid->v[i][0] += dt * Fluid->a[i][0];
////			Fluid->v[i][1] += dt * Fluid->a[i][1];
////			Fluid->v[i][2] += dt * Fluid->a[i][2];
////
////			Fluid->x[i][0] += dt * Fluid->v[i][0];
////			Fluid->x[i][1] += dt * Fluid->v[i][1];
////			Fluid->x[i][2] += dt * Fluid->v[i][2];
////		};
////		auto end = std::chrono::high_resolution_clock::now();
////		duration3 = end - start;
////		std::cout << "for Time: " << duration3.count() << " ms`" << std::endl;
////	}
//
//	{
////		auto start = std::chrono::high_resolution_clock::now();
////		size_t patch = size / 16;
//		size_t left = size % 16;
//
//		for (int i = 0; i < size - left; i += 16)
//		{
//			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->a[i][0], dt, 16 * 3);
//			VM_Math::addSIMD(&Fluid->v[i][0], &Fluid->v[i][0], &Fluid->tmpv[i][0], 16 * 3);
//			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->v[i][0], dt, 16 * 3);
//			VM_Math::addSIMD(&Fluid->x[i][0], &Fluid->x[i][0], &Fluid->tmpv[i][0], 16 * 3);
//		};
//
//		for (size_t i = size - left; i < size; ++i)
//		{
//			Fluid->v[i][0] += dt * Fluid->a[i][0];
//			Fluid->v[i][1] += dt * Fluid->a[i][1];
//			Fluid->v[i][2] += dt * Fluid->a[i][2];
//
//			Fluid->x[i][0] += dt * Fluid->v[i][0];
//			Fluid->x[i][1] += dt * Fluid->v[i][1];
//			Fluid->x[i][2] += dt * Fluid->v[i][2];
//		}
//
////		auto end = std::chrono::high_resolution_clock::now();
////		duration4 = end - start;
////		std::cout << "SIMD Time: " << duration4.count() << " ms" << std::endl;
//	}
//
////	std::cout << "Acc Rate 12: " << duration1.count() / duration2.count() << std::endl;
////	std::cout << "Acc Rate 13: " << duration1.count() / duration3.count() << std::endl;
////	std::cout << "Acc Rate 23: " << duration2.count() / duration3.count() << std::endl;
////	std::cout << "Acc Rate 14: " << duration1.count() / duration4.count() << std::endl;
//}

#ifdef TEST_DFSPH
#include <iostream>
#include <random>
#include <GEO/GEO_BVH.h>
int main()
{
	int N = 2 << 20;
	std::vector<float> a(N);
	std::vector<float> x(N);
	std::vector<float> b(N);
	// random init
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0, 1);
	for (int i = 0; i < N; ++i)
	{
		a[i] = dis(gen);
		x[i] = dis(gen);
		b[i] = dis(gen);
	}
	for (int i = 0; i < N; i += 16)
	{
		VM_Math::mulSIMD(&b[i], &a[i], 0.1f, 16);
		VM_Math::addSIMD(&x[i], &x[i], &b[i], 16);
	}
	std::cout << "size: " << N << std::endl;
	std::cout << x[0] << std::endl;
	return 0;
}
#endif
