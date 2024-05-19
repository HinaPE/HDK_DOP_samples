#include "DFSPH.h"
#include <GU/GU_Detail.h>
#include <GU/GU_NeighbourList.h>
#include <UT/UT_ParallelUtil.h>

using Vector3 = HinaPE::TBB::float3;

struct Cubic
{
	inline static void set_radius(float r)
	{
		_r = r;
		const float pi = 3.14159265358979323846f;
		const float h3 = _r * _r * _r;
		_k = 8.f / (pi * h3);
		_l = 48.f / (pi * h3);
		_W_0 = W(Vector3{0, 0, 0});
	}
	inline static float W(const float r)
	{
		float res = 0.0;
		const float q = r / _r;
		if (q <= 1.0)
		{
			if (q <= 0.5)
			{
				const float q2 = q * q;
				const float q3 = q2 * q;
				res = _k * (6.f * q3 - 6.f * q2 + 1.f);
			} else
			{
				res = _k * (2.f * pow(1.f - q, 3.f));
			}
		}
		return res;
	}
	inline static float W(const Vector3 &r)
	{
		return W(r.length());
	}

	inline static Vector3 gradW(const Vector3 &r)
	{
		Vector3 res;
		const float rl = r.length();
		const float q = rl / _r;
		if ((rl > 1.0e-9) && (q <= 1.0))
		{
			Vector3 gradq = r / rl;
			gradq /= _r;
			if (q <= 0.5)
			{
				res = _l * q * (3.f * q - 2.f) * gradq;
			} else
			{
				const float factor = 1.f - q;
				res = _l * (-factor * factor) * gradq;
			}
		} else
			res = Vector3{0, 0, 0};

		return res;
	}

	inline static float W_zero()
	{
		return _W_0;
	}

	inline static float _r;
	inline static float _k;
	inline static float _l;
	inline static float _W_0;
};

using Kernel = Cubic;

inline void parallel_for(size_t n, const std::function<void(size_t)> &f)
{
	UTparallelForEachNumber((int) n, [&](const UT_BlockedRange<int> &range)
	{
		for (size_t i = range.begin(); i != range.end(); ++i) { f(i); }
	});
}

HinaPE::TBB::DFSPH::DFSPH(float _kernel_radius) : size(0), kr(_kernel_radius)
{
	Fluid = std::make_shared<FluidData>();
	Kernel::set_radius(_kernel_radius);
}
void HinaPE::TBB::DFSPH::solve(float dt, const GU_Detail *gdp)
{
	printf("dt: %f\n", dt);
	// ==================== 1. Refresh Fluid Data ====================
	static std::vector<UT_Array<GA_Offset>> NL;
	if (size < gdp->getNumPoints())
	{
		auto new_size = gdp->getNumPoints();
		Fluid->x.resize(new_size, {0, 0, 0});
		Fluid->v.resize(new_size, {0, 0, 0});
		Fluid->a.resize(new_size, {0, 0, 0});
		Fluid->m.resize(new_size, DEFAULT_M);
		Fluid->V.resize(new_size, DEFAULT_V);
		Fluid->rho.resize(new_size, 0);
		Fluid->factor.resize(new_size, 0);
		Fluid->density_adv.resize(new_size, 0);
		Fluid->nn.resize(new_size, 0);
		Fluid->tmp.resize(new_size, 0);
		size = new_size;

		{
			GA_Offset pt_off;
			GA_FOR_ALL_PTOFF(gdp, pt_off)
				{
					GA_Index pt_idx = gdp->pointIndex(pt_off);
					UT_Vector3 pos = gdp->getPos3(pt_off);
					Fluid->x[pt_idx] = float3{pos.x(), pos.y(), pos.z()};
				}
		}

		NL.resize(new_size);
	}



	// ==================== 2. Build Neighbors ====================
	GU_NeighbourList Searcher;
	GU_NeighbourListParms NLP;
	NLP.setRadius(kr);
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
	parallel_for(size, [&](size_t i)
	{
		Fluid->rho[i] = Fluid->V[i] * Kernel::W_zero();
		for (auto j: NL[i])
			Fluid->rho[i] += Fluid->V[j] * Kernel::W(Fluid->x[i] - Fluid->x[j]);
		Fluid->rho[i] *= REST_DENSITY;

		float sum_grad_p_k = 0.;
		float3 grad_p_i{0, 0, 0};
		for (auto j: NL[i])
		{
			float3 grad_p_j = -Fluid->V[j] * Kernel::gradW(Fluid->x[i] - Fluid->x[j]);
			sum_grad_p_k += grad_p_j.length2();
			grad_p_i -= grad_p_j;
		}
		sum_grad_p_k += grad_p_i.length2();
		if (sum_grad_p_k > 1e-6)
			Fluid->factor[i] = -1.f / sum_grad_p_k;
		else
			Fluid->factor[i] = 0;
	});



	// ==================== 4. Divergence Solver ====================
	parallel_for(size, [&](size_t i) // compute density change
	{
		float density_adv = 0;
		for (auto j: NL[i])
			density_adv += Fluid->V[j] * (Fluid->v[i] - Fluid->v[j]).dot(Kernel::gradW(Fluid->x[i] - Fluid->x[j]));
		Fluid->density_adv[i] = std::max(density_adv, 0.f);
	});
	size_t iter = 0;
	float avg_density_err = 0.0;
	while (iter < 1 || iter < MAX_ITERATIONS)
	{
		parallel_for(size, [&](size_t i)
		{
			float b_i = Fluid->density_adv[i];
			float k_i = b_i * Fluid->factor[i] / dt;
			float3 dv{0, 0, 0};
			for (auto j: NL[i])
			{
				float b_j = Fluid->density_adv[j];
				float k_j = b_j * Fluid->factor[j] / dt;
				float k_sum = k_i + k_j; // for multiphase fluid
				if (std::abs(k_sum) > 1e-5)
				{
					float3 grad_p_j = -Fluid->V[j] * Kernel::gradW(Fluid->x[i] - Fluid->x[j]);
					dv -= dt * k_sum * grad_p_j;
				}
			}
			Fluid->v[i] += dv;
		});
		parallel_for(size, [&](size_t i) // compute density change
		{
			float density_adv = 0;
			for (auto j: NL[i])
				density_adv += Fluid->V[j] * (Fluid->v[i] - Fluid->v[j]).dot(Kernel::gradW(Fluid->x[i] - Fluid->x[j]));
			Fluid->density_adv[i] = std::max(density_adv, 0.f);
		});
		avg_density_err = 0;
		for (size_t i = 0; i < size; ++i)
			avg_density_err += REST_DENSITY * Fluid->density_adv[i];
		avg_density_err /= size;

		float eta = 1.f / dt * .1f * 0.01f * REST_DENSITY;
		if (avg_density_err <= eta)
			break;
		++iter;
	}
	printf("(div solver) iter: %zu, avg_density_err: %f\n", iter, avg_density_err);



	// ==================== 5. Non-Pressure Force and Predict Velocity ====================
	constexpr float d = 10;
	const float diameter = 2 * PARTICLE_RADIUS;
	const float diameter2 = diameter * diameter;
	const float kr2 = kr * kr;
	parallel_for(size, [&](size_t i)
	{
		float3 dv = GRAVITY;
		for (auto j: NL[i])
		{
			// Surface Tension
			const float3 r = Fluid->x[i] - Fluid->x[j];
			const float r2 = r.dot(r);
			const float r1 = std::sqrt(r2);
			if (r2 > diameter2)
				dv -= SURFACE_TENSION / Fluid->m[i] * Fluid->m[j] * r * Kernel::W(r1);
			else
				dv -= SURFACE_TENSION / Fluid->m[i] * Fluid->m[j] * r * Kernel::W(diameter);

			// Fluid Viscosity
			float v_xy = (Fluid->v[i] - Fluid->v[j]).dot(r);
			float3 f_v = d * VISCOSITY * (Fluid->m[j] / (Fluid->rho[j])) * v_xy / (r2 + 0.01f * kr2) * Kernel::gradW(r);
			dv += f_v;
		}
		Fluid->a[i] = dv;
		Fluid->v[i] += dt * Fluid->a[i];
	});



	// ==================== 6. Pressure Solver ====================
	parallel_for(size, [&](size_t i) // compute density adv
	{
		float delta = 0;
		for (auto j: NL[i])
			delta += Fluid->V[j] * (Fluid->v[i] - Fluid->v[j]).dot(Kernel::gradW(Fluid->x[i] - Fluid->x[j]));
		float density_adv = Fluid->rho[i] / REST_DENSITY + dt * delta;
		Fluid->density_adv[i] = std::max(density_adv, 1.f);
	});
	iter = 0;
	avg_density_err = 0.0;
	while (iter < 1 || iter < MAX_ITERATIONS)
	{
		parallel_for(size, [&](size_t i)
		{
			float b_i = Fluid->density_adv[i] - 1.f;
			float k_i = b_i * Fluid->factor[i] / (dt * dt);
			float3 dv{0, 0, 0};
			for (auto j: NL[i])
			{
				float b_j = Fluid->density_adv[j] - 1.f;
				float k_j = b_j * Fluid->factor[j] / (dt * dt);
				float k_sum = k_i + k_j;
				if (std::abs(k_sum) > 1e-5)
				{
					float3 grad_p_j = -Fluid->V[j] * Kernel::gradW(Fluid->x[i] - Fluid->x[j]);
					dv -= dt * k_sum * grad_p_j;
				}
			}
			Fluid->v[i] += dv;
		});
		parallel_for(size, [&](size_t i) // compute density adv
		{
			float delta = 0;
			for (auto j: NL[i])
				delta += Fluid->V[j] * (Fluid->v[i] - Fluid->v[j]).dot(Kernel::gradW(Fluid->x[i] - Fluid->x[j]));
			float density_adv = Fluid->rho[i] / REST_DENSITY + dt * delta;
			Fluid->density_adv[i] = std::max(density_adv, 1.f);
		});
		avg_density_err = 0;
		for (size_t i = 0; i < size; ++i)
			avg_density_err += REST_DENSITY * (Fluid->density_adv[i] - 1.f);
		avg_density_err /= size;
		float eta = 0.05f * 0.01f * REST_DENSITY;
		if (avg_density_err <= eta)
			break;
		++iter;
	}
	printf("(pressure solver) iter: %zu, avg_density_err: %f\n", iter, avg_density_err);



	// ==================== 7. Advect ====================
	parallel_for(size, [&](size_t i)
	{
		Fluid->x[i] += dt * Fluid->v[i];
	});


	// ==================== 8. Enforce Boundary ====================
	parallel_for(size, [&](size_t i)
	{
		float3 collision_normal{0, 0, 0};
		if (Fluid->x[i].x() > MaxBound.x())
		{
			Fluid->x[i].x() = MaxBound.x();
			collision_normal.x() += 1;
		}
		if (Fluid->x[i].x() < -MaxBound.x())
		{
			Fluid->x[i].x() = -MaxBound.x();
			collision_normal.x() -= 1;
		}
		if (!TOP_OPEN)
		{
			if (Fluid->x[i].y() > MaxBound.y())
			{
				Fluid->x[i].y() = MaxBound.y();
				collision_normal.y() += 1;
			}
		}
		if (Fluid->x[i].y() < -MaxBound.y())
		{
			Fluid->x[i].y() = -MaxBound.y();
			collision_normal.y() -= 1;
		}
		if (Fluid->x[i].z() > MaxBound.z())
		{
			Fluid->x[i].z() = MaxBound.z();
			collision_normal.z() += 1;
		}
		if (Fluid->x[i].z() < -MaxBound.z())
		{
			Fluid->x[i].z() = -MaxBound.z();
			collision_normal.z() -= 1;
		}
		collision_normal.normalize();
		constexpr float c_f = 0.5f;
		Fluid->v[i] -= (1.f + c_f) * Fluid->v[i].dot(collision_normal) * collision_normal;
	});
}

#ifdef TEST_DFSPH
#include <iostream>
int main()
{
	std::cout << "Hello, TBB!" << std::endl;
	return 0;
}
#endif
