#include "DFSPH.h"

#include <UT/UT_ParallelUtil.h>
#include <GU/GU_NeighbourList.h>
#include <VM/VM_Math.h>
#include <iostream>
#include <chrono>

inline void parallel_for(size_t n, const std::function<void(size_t)> &f)
{
	UTparallelForEachNumber((int) n, [&](const UT_BlockedRange<int> &range)
	{
		for (size_t i = range.begin(); i != range.end(); ++i) { f(i); }
	});
}

HinaPE::SIMD::DFSPH::DFSPH(float _kernel_radius)
{
	Fluid = std::make_shared<FluidSIMD>();

	Searcher = std::make_shared<GU_NeighbourList>();
	nlp = std::make_shared<GU_NeighbourListParms>();
	nlp->setRadius(_kernel_radius);
	nlp->setOverrideRadius(true);
	nlp->setMode(GU_NeighbourListParms::InteractionMode::UNIFORM);
}
void HinaPE::SIMD::DFSPH::resize(size_t n)
{
	if (size == n)
		return;
	Fluid->x.resize(n, {0, 0, 0});
	Fluid->v.resize(n, {0, 0, 0});
	Fluid->a.resize(n, {0, 0, 0});
	Fluid->m.resize(n, 1000.f * 0.9f * 0.02f * 0.02f * 0.02f);
	Fluid->V.resize(n, 0.9f * 0.02f * 0.02f * 0.02f);
	Fluid->rho.resize(n, 0);
	Fluid->factor.resize(n, 0);
	Fluid->density_adv.resize(n, 0);
	Fluid->nn.resize(n, 0);
	Fluid->tmpv.resize(n, {0, 0, 0});
	Fluid->tmp.resize(n, 0);
	size = n;
}
void HinaPE::SIMD::DFSPH::solve(float dt, GU_Detail &gdp)
{
	// ==================== 1. Build Neighbors ====================
	Searcher->build(&gdp, *nlp);



	// ==================== 2. Compute Density and Factor ====================
	parallel_for(size, [&](size_t i)
	{
		UT_Array<GA_Offset> neighbor_list;
		Searcher->getNeighbours(i, &gdp, neighbor_list);
		Fluid->nn[i] = neighbor_list.size();
	});

	std::transform(Fluid->a.begin(), Fluid->a.end(), Fluid->a.begin(), [](std::array<float, 3> a) { return std::array<float, 3>{0, -9.8, 0}; });

//	std::chrono::duration<double, std::milli> duration1, duration2, duration3, duration4;
//	{
//		auto start = std::chrono::high_resolution_clock::now();
//		size_t patch = size / 16;
//		size_t left = size % 16;
//
//		parallel_for(patch, [&](size_t p)
//		{
//			size_t i = p * 16;
//			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->a[i][0], dt, 16 * 3);
//			VM_Math::addSIMD(&Fluid->v[i][0], &Fluid->v[i][0], &Fluid->tmpv[i][0], 16 * 3);
//			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->v[i][0], dt, 16 * 3);
//			VM_Math::addSIMD(&Fluid->x[i][0], &Fluid->x[i][0], &Fluid->tmpv[i][0], 16 * 3);
//		});
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
//		auto end = std::chrono::high_resolution_clock::now();
//		duration1 = end - start;
//		std::cout << "parallel + SIMD Time: " << duration1.count() << " ms" << std::endl;
//	}
//
//	{
//		auto start = std::chrono::high_resolution_clock::now();
//		parallel_for(size, [&](size_t i)
//		{
//			Fluid->v[i][0] += dt * Fluid->a[i][0];
//			Fluid->v[i][1] += dt * Fluid->a[i][1];
//			Fluid->v[i][2] += dt * Fluid->a[i][2];
//
//			Fluid->x[i][0] += dt * Fluid->v[i][0];
//			Fluid->x[i][1] += dt * Fluid->v[i][1];
//			Fluid->x[i][2] += dt * Fluid->v[i][2];
//		});
//		auto end = std::chrono::high_resolution_clock::now();
//		duration2 = end - start;
//		std::cout << "parallel Time: " << duration2.count() << " ms" << std::endl;
//	}
//
//	{
//		auto start = std::chrono::high_resolution_clock::now();
//		for (int i = 0; i < size; ++i)
//		{
//			Fluid->v[i][0] += dt * Fluid->a[i][0];
//			Fluid->v[i][1] += dt * Fluid->a[i][1];
//			Fluid->v[i][2] += dt * Fluid->a[i][2];
//
//			Fluid->x[i][0] += dt * Fluid->v[i][0];
//			Fluid->x[i][1] += dt * Fluid->v[i][1];
//			Fluid->x[i][2] += dt * Fluid->v[i][2];
//		};
//		auto end = std::chrono::high_resolution_clock::now();
//		duration3 = end - start;
//		std::cout << "for Time: " << duration3.count() << " ms`" << std::endl;
//	}

	{
//		auto start = std::chrono::high_resolution_clock::now();
//		size_t patch = size / 16;
		size_t left = size % 16;

		for (int i = 0; i < size - left; i += 16)
		{
			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->a[i][0], dt, 16 * 3);
			VM_Math::addSIMD(&Fluid->v[i][0], &Fluid->v[i][0], &Fluid->tmpv[i][0], 16 * 3);
			VM_Math::mulSIMD(&Fluid->tmpv[i][0], &Fluid->v[i][0], dt, 16 * 3);
			VM_Math::addSIMD(&Fluid->x[i][0], &Fluid->x[i][0], &Fluid->tmpv[i][0], 16 * 3);
		};

		for (size_t i = size - left; i < size; ++i)
		{
			Fluid->v[i][0] += dt * Fluid->a[i][0];
			Fluid->v[i][1] += dt * Fluid->a[i][1];
			Fluid->v[i][2] += dt * Fluid->a[i][2];

			Fluid->x[i][0] += dt * Fluid->v[i][0];
			Fluid->x[i][1] += dt * Fluid->v[i][1];
			Fluid->x[i][2] += dt * Fluid->v[i][2];
		}

//		auto end = std::chrono::high_resolution_clock::now();
//		duration4 = end - start;
//		std::cout << "SIMD Time: " << duration4.count() << " ms" << std::endl;
	}

//	std::cout << "Acc Rate 12: " << duration1.count() / duration2.count() << std::endl;
//	std::cout << "Acc Rate 13: " << duration1.count() / duration3.count() << std::endl;
//	std::cout << "Acc Rate 23: " << duration2.count() / duration3.count() << std::endl;
//	std::cout << "Acc Rate 14: " << duration1.count() / duration4.count() << std::endl;
}

#ifdef TEST_DFSPH
#include <iostream>
int main()
{
	return 0;
}
#endif
