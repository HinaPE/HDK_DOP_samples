#include "DFSPH.cuh"
#include <cuda_runtime.h>
#include "common/helper_cuda.h"
#include "common/helper_math.h"
#include "neighbor/include/cuNSearch.h"
#include "neighbor/src/cuNSearchDeviceData.h"
#include <thrust/execution_policy.h>

__constant__ float GPU_PARTICLE_RADIUS;
__constant__ float GPU_KERNEL_RADIUS;
__constant__ float GPU_KERNEL_K;
__constant__ float GPU_KERNEL_L;
__constant__ float GPU_REST_DENSITY;
__constant__ float GPU_VISCOSITY;
__constant__ float GPU_SURFACE_TENSION;
__constant__ float3 GPU_MAX_BOUND;
__constant__ bool GPU_TOP_OPEN;
static inline __device__ float W(const float r)
{
	float res = 0.f;
	const float q = r / GPU_KERNEL_RADIUS;
	if (q <= 1.f)
	{
		if (q <= 0.5f)
		{
			const float q2 = q * q;
			const float q3 = q2 * q;
			res = GPU_KERNEL_K * (6.f * q3 - 6.f * q2 + 1.f);
		} else
		{
			res = GPU_KERNEL_K * (2.f * powf(1.f - q, 3.f));
		}
	}
	return res;
}
static inline __device__ float W(const float3 r) { return W(length(r)); }
static inline __device__ float3 gradW(const float3 r)
{
	float3 res;
	const float rl = length(r);
	const float q = rl / GPU_KERNEL_RADIUS;
	if ((rl > 1.0e-9) && (q <= 1.f))
	{
		float3 gradq = {r.x / rl, r.y / rl, r.z / rl};
		gradq.x /= GPU_KERNEL_RADIUS;
		gradq.y /= GPU_KERNEL_RADIUS;
		gradq.z /= GPU_KERNEL_RADIUS;
		if (q <= 0.5f)
		{
			res = GPU_KERNEL_L * q * (3.f * q - 2.f) * gradq;
		} else
		{
			const float factor = 1.f - q;
			res = GPU_KERNEL_L * (-factor * factor) * gradq;
		}
	} else
		res = {0.f, 0.f, 0.f};
	return res;
}
static inline __device__ float W_zero() { return W(make_float3(0.f, 0.f, 0.f)); }

HinaPE::CUDA::DFSPH::DFSPH() : size(0)
{
	Fluid = std::make_shared<FluidGPU>();
	Searcher = std::make_shared<cuNSearch::NeighborhoodSearch>(KERNEL_RADIUS);
	fluid_idx = Searcher->add_point_set(&(Fluid->x.data()->x), Fluid->x.size(), true, true, true);
	Searcher->set_active(true);

	cudaMalloc((void **) &GPU_PARTICLE_RADIUS, sizeof(float));
	cudaMalloc((void **) &GPU_KERNEL_RADIUS, sizeof(float));
	cudaMalloc((void **) &GPU_KERNEL_K, sizeof(float));
	cudaMalloc((void **) &GPU_KERNEL_L, sizeof(float));
	cudaMalloc((void **) &GPU_REST_DENSITY, sizeof(float));
	cudaMalloc((void **) &GPU_VISCOSITY, sizeof(float));
	cudaMalloc((void **) &GPU_SURFACE_TENSION, sizeof(float));
	cudaMalloc((void **) &GPU_MAX_BOUND, sizeof(float3));
	cudaMalloc((void **) &GPU_TOP_OPEN, sizeof(bool));

	set_gpu_constants();
}

void HinaPE::CUDA::DFSPH::set_gpu_constants()
{
	float k = 8.f / (3.14159265358979323846f * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS);
	float l = 48.f / (3.14159265358979323846f * KERNEL_RADIUS * KERNEL_RADIUS * KERNEL_RADIUS);
	cudaMemcpyToSymbol(GPU_PARTICLE_RADIUS, &PARTICLE_RADIUS, sizeof(float));
	cudaMemcpyToSymbol(GPU_KERNEL_RADIUS, &KERNEL_RADIUS, sizeof(float));
	cudaMemcpyToSymbol(GPU_KERNEL_K, &k, sizeof(float));
	cudaMemcpyToSymbol(GPU_KERNEL_L, &l, sizeof(float));
	cudaMemcpyToSymbol(GPU_REST_DENSITY, &REST_DENSITY, sizeof(float));
	cudaMemcpyToSymbol(GPU_SURFACE_TENSION, &SURFACE_TENSION, sizeof(float));
	cudaMemcpyToSymbol(GPU_VISCOSITY, &VISCOSITY, sizeof(float));
	cudaMemcpyToSymbol(GPU_MAX_BOUND, &MaxBound, 3 * sizeof(float));
	cudaMemcpyToSymbol(GPU_TOP_OPEN, &TOP_OPEN, sizeof(bool));
}

void HinaPE::CUDA::DFSPH::resize(size_t n)
{
	// ==================== 1. Refresh Fluid Data ====================
	if (size < n)
	{
		auto new_size = n;
		Fluid->x.resize(new_size, {0, 0, 0}); // NO NEED TO RESIZE
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

		Searcher->resize_point_set(fluid_idx, &(Fluid->x.data()->x), size);
	}
}

void HinaPE::CUDA::DFSPH::solve(float dt)
{
	// ==================== 2. Build Neighbors ====================
	Searcher->update_point_set(fluid_idx);
	Searcher->find_neighbors();
	cuNSearch::PointSet::NeighborSet &neighbor_set = Searcher->point_set(fluid_idx).get_raw_neighbor_set(fluid_idx);
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					nn = Fluid->nn.data(),
					dNeighbors = neighbor_set.d_Neighbors.data(),
					dCounts = neighbor_set.d_NeighborCounts.data(),
					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
			] __device__(size_t i)
			{
				nn[i] = dCounts[i];
			});



	// ==================== 3. Compute Density and Factor ====================
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					V = Fluid->V.data(),
					rho = Fluid->rho.data(),
					factor = Fluid->factor.data(),
					nn = Fluid->nn.data(),
					dNeighbors = neighbor_set.d_Neighbors.data(),
					dCounts = neighbor_set.d_NeighborCounts.data(),
					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
			] __device__(size_t i)
			{
				float sum_grad = 0.;
				float3 grad_i{0, 0, 0};
				rho[i] = V[i] * W_zero();
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];
					rho[i] += V[j] * W(x[i] - x[j]);
					float3 grad_j = -V[j] * gradW(x[i] - x[j]);
					sum_grad += dot(grad_j, grad_j);
					grad_i -= grad_j;
				}
				rho[i] *= GPU_REST_DENSITY;
				sum_grad += dot(grad_i, grad_i);
				if (sum_grad > 1e-6)
					factor[i] = -1.f / sum_grad;
				else
					factor[i] = 0;
			});



	// ==================== 4. Divergence Solver ====================
	thrust::for_each( // compute density change
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data(),
					V = Fluid->V.data(),
					density_adv = Fluid->density_adv.data(),
					dNeighbors = neighbor_set.d_Neighbors.data(),
					dCounts = neighbor_set.d_NeighborCounts.data(),
					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
			] __device__(size_t i)
			{
				float _density_adv = 0;
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];
					_density_adv += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
				}
				density_adv[i] = std::max(_density_adv, 0.f);
			});
	DIVERGENCE_ITERS = 0;
	float avg_density_err = 0.0;
	while (DIVERGENCE_ITERS < 1 || DIVERGENCE_ITERS < MAX_ITERATIONS)
	{
		thrust::for_each(
				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
				[
						x = Fluid->x.data(),
						v = Fluid->v.data(),
						V = Fluid->V.data(),
						factor = Fluid->factor.data(),
						density_adv = Fluid->density_adv.data(),
						dNeighbors = neighbor_set.d_Neighbors.data(),
						dCounts = neighbor_set.d_NeighborCounts.data(),
						dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
						dt
				] __device__(size_t i)
				{
					float b_i = density_adv[i];
					float k_i = b_i * factor[i] / dt;
					float3 dv{0, 0, 0};
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						float b_j = density_adv[j];
						float k_j = b_j * factor[j] / dt;
						float k_sum = k_i + k_j;
						if (std::abs(k_sum) > 1e-5)
						{
							float3 grad_j = -V[j] * gradW(x[i] - x[j]);
							dv -= dt * k_sum * grad_j;
						}
					}
					v[i] += dv;
				});
		thrust::for_each( // compute density change
				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
				[
						x = Fluid->x.data(),
						v = Fluid->v.data(),
						V = Fluid->V.data(),
						density_adv = Fluid->density_adv.data(),
						dNeighbors = neighbor_set.d_Neighbors.data(),
						dCounts = neighbor_set.d_NeighborCounts.data(),
						dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
				] __device__(size_t i)
				{
					float _density_adv = 0;
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						_density_adv += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
					}
					density_adv[i] = std::max(_density_adv, 0.f);
				});
		avg_density_err = 0;
		for (size_t i = 0; i < size; ++i)
			avg_density_err += REST_DENSITY * Fluid->density_adv[i];
		avg_density_err /= size;

		float eta = 1.f / dt * .1f * 0.01f * REST_DENSITY;
		if (avg_density_err <= eta)
			break;
		++DIVERGENCE_ITERS;
	}



	// ==================== 5. Non-Pressure Force and Predict Velocity ====================
	constexpr float d = 10;
	const float diameter = 2 * PARTICLE_RADIUS;
	const float diameter2 = diameter * diameter;
	const float kr2 = KERNEL_RADIUS * KERNEL_RADIUS;
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data(),
					V = Fluid->V.data(),
					m = Fluid->m.data(),
					a = Fluid->a.data(),
					rho = Fluid->rho.data(),
					density_adv = Fluid->density_adv.data(),
					dNeighbors = neighbor_set.d_Neighbors.data(),
					dCounts = neighbor_set.d_NeighborCounts.data(),
					dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
					diameter,
					diameter2,
					kr2
			] __device__(size_t i)
			{
				float3 dv{0, -9.8f, 0};
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];

					// Surface Tension
					const float3 r = x[i] - x[j];
					const float r2 = dot(r, r);
					const float r1 = std::sqrt(r2);
					if (r2 > diameter2)
						dv -= GPU_SURFACE_TENSION / m[i] * m[j] * r * W(r1);
					else
						dv -= GPU_SURFACE_TENSION / m[i] * m[j] * r * W(diameter);

					// Fluid Viscosity
					float v_xy = dot(v[i] - v[j], r);
					float3 f_v = d * GPU_VISCOSITY * (m[j] / (rho[j])) * v_xy / (r2 + 0.01f * kr2) * gradW(r);
					dv += f_v;
				}
				a[i] = dv;
			});
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					v = Fluid->v.data(),
					a = Fluid->a.data(),
					dt
			] __device__(size_t i)
			{
				v[i] += dt * a[i];
			});



	// ==================== 6. Pressure Solver ====================
	thrust::for_each( // compute density adv
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data(),
					V = Fluid->V.data(),
					rho = Fluid->rho.data(),
					density_adv = Fluid->density_adv.data(),
					dNeighbors = neighbor_set.d_Neighbors.data(),
					dCounts = neighbor_set.d_NeighborCounts.data(),
					dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
					dt
			] __device__(size_t i)
	{
		float delta = 0;
		for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
		{
			uint j = dNeighbors[dOffsets[i] + _n_idx];
			delta += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
		}
		float _density_adv = rho[i] / GPU_REST_DENSITY + dt * delta;
		density_adv[i] = std::max(_density_adv, 1.f);
	});
	PRESSURE_ITERS = 0;
	avg_density_err = 0.0;
	while (PRESSURE_ITERS < 1 || PRESSURE_ITERS < MAX_ITERATIONS)
	{
		thrust::for_each(
				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
				[
						x = Fluid->x.data(),
						v = Fluid->v.data(),
						V = Fluid->V.data(),
						factor = Fluid->factor.data(),
						density_adv = Fluid->density_adv.data(),
						dNeighbors = neighbor_set.d_Neighbors.data(),
						dCounts = neighbor_set.d_NeighborCounts.data(),
						dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
						dt
				] __device__(size_t i)
		{
			float b_i = density_adv[i] - 1.f;
			float k_i = b_i * factor[i] / (dt * dt);
			float3 dv{0, 0, 0};
			for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
			{
				uint j = dNeighbors[dOffsets[i] + _n_idx];
				float b_j = density_adv[j] - 1.f;
				float k_j = b_j * factor[j] / (dt * dt);
				float k_sum = k_i + k_j;
				if (std::abs(k_sum) > 1e-5)
				{
					float3 grad_p_j = -V[j] * gradW(x[i] - x[j]);
					dv -= dt * k_sum * grad_p_j;
				}
			}
			v[i] += dv;
		});
		thrust::for_each( // compute density adv
				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
				[
						x = Fluid->x.data(),
						v = Fluid->v.data(),
						V = Fluid->V.data(),
						rho = Fluid->rho.data(),
						density_adv = Fluid->density_adv.data(),
						dNeighbors = neighbor_set.d_Neighbors.data(),
						dCounts = neighbor_set.d_NeighborCounts.data(),
						dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
						dt
				] __device__(size_t i)
				{
					float delta = 0;
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						delta += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
					}
					float _density_adv = rho[i] / GPU_REST_DENSITY + dt * delta;
					density_adv[i] = std::max(_density_adv, 1.f);
				});
		avg_density_err = 0;
		for (size_t i = 0; i < size; ++i)
			avg_density_err += REST_DENSITY * (Fluid->density_adv[i] - 1.f);
		avg_density_err /= size;

		float eta = 0.05f * 0.01f * REST_DENSITY;
		if (avg_density_err <= eta)
			break;
		++PRESSURE_ITERS;
	}



	// ==================== 7. Advect ====================
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data(),
					dt
			] __device__(size_t i)
	{
		x[i] += dt * v[i];
	});



	// ==================== 8. Enforce Boundary ====================
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data()
			] __device__(size_t i)
	{
		float3 normal{0, 0, 0};
		if (x[i].x > GPU_MAX_BOUND.x)
		{
			x[i].x = GPU_MAX_BOUND.x;
			normal.x += 1;
		}
		if (x[i].x < -GPU_MAX_BOUND.x)
		{
			x[i].x = -GPU_MAX_BOUND.x;
			normal.x -= 1;
		}
		if (!GPU_TOP_OPEN)
		{
			if (x[i].y > GPU_MAX_BOUND.y)
			{
				x[i].y = GPU_MAX_BOUND.y;
				normal.y += 1;
			}
		}
		if (x[i].y < -GPU_MAX_BOUND.y)
		{
			x[i].y = -GPU_MAX_BOUND.y;
			normal.y -= 1;
		}
		if (x[i].z > GPU_MAX_BOUND.z)
		{
			x[i].z = GPU_MAX_BOUND.z;
			normal.z += 1;
		}
		if (x[i].z < -GPU_MAX_BOUND.z)
		{
			x[i].z = -GPU_MAX_BOUND.z;
			normal.z -= 1;
		}
		if (length(normal) > std::numeric_limits<float>::epsilon())
		{
			normal = normalize(normal);
			constexpr float c_f = 0.5f;
			v[i] -= (1.f + c_f) * dot(v[i], normal) * normal;
		}
	});
}

//void HinaPE::CUDA::DFSPH::solve(float dt)
//{
//
//	// ==================== 1. Build Neighbors ====================
//	if (need_reload)
//	{
//		Searcher->resize_point_set(fluid_idx, &(Fluid->x.data()->x), size);
//		need_reload = false;
//	}
//	Searcher->update_point_set(fluid_idx);
//	Searcher->find_neighbors();
//	cuNSearch::PointSet::NeighborSet &neighbor_set = Searcher->point_set(fluid_idx).get_raw_neighbor_set(fluid_idx);
//
//
//
//	// ==================== 2. Compute Density and Factor ====================
//	thrust::for_each(
//			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//			[
//					x = Fluid->x.data(),
//					V = Fluid->V.data(),
//					rho = Fluid->rho.data(),
//					factor = Fluid->factor.data(),
//					nn = Fluid->nn.data(),
//					dNeighbors = neighbor_set.d_Neighbors.data(),
//					dCounts = neighbor_set.d_NeighborCounts.data(),
//					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
//			] __device__(size_t i)
//			{
//				float _rho_i = V[i] * W0();
//				float _sum_grad_p_k = 0;
//				float3 _grad_p_i = {0.f, 0.f, 0.f};
//				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//				{
//					uint j = dNeighbors[dOffsets[i] + _n_idx];
//					_rho_i += V[j] * W(x[i] - x[j]);
//
//					float3 _grad_p_j = -V[j] * gradW(x[i] - x[j]);
//					_sum_grad_p_k += dot(_grad_p_j, _grad_p_j);
//					_grad_p_i -= _grad_p_j;
//				}
//
//				_sum_grad_p_k += dot(_grad_p_i, _grad_p_i);
//				rho[i] = _rho_i * REST_DENSITY;
//				if (_sum_grad_p_k > 1e-6f)
//					factor[i] = -1.f / _sum_grad_p_k;
//				else
//					factor[i] = 0;
//
//				nn[i] = dCounts[i];
//			});
//
//
//
//	// ==================== 3. Divergence Solve ====================
//	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
//	_) { return _ / dt; });
//	thrust::for_each(
//			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//			[
//					x = Fluid->x.data(),
//					v = Fluid->v.data(),
//					V = Fluid->V.data(),
//					density_adv = Fluid->density_adv.data(),
//					dNeighbors = neighbor_set.d_Neighbors.data(),
//					dCounts = neighbor_set.d_NeighborCounts.data(),
//					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
//			] __device__(size_t i)
//			{
//				float _d_adv = 0;
//				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//				{
//					uint j = dNeighbors[dOffsets[i] + _n_idx];
//					_d_adv += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
//				}
//				density_adv[i] = max(_d_adv, 0.f);
//			});
//	uint iteration_v = 0;
//	float avg_density_error = 0;
//	while (iteration_v < 1 || iteration_v < 100)
//	{
//		thrust::for_each(
//				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//				[
//						x = Fluid->x.data(),
//						v = Fluid->v.data(),
//						V = Fluid->V.data(),
//						factor = Fluid->factor.data(),
//						density_adv = Fluid->density_adv.data(),
//						dNeighbors = neighbor_set.d_Neighbors.data(),
//						dCounts = neighbor_set.d_NeighborCounts.data(),
//						dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
//						dt
//				] __device__(size_t i)
//				{
//					float _k_i = density_adv[i] * factor[i];
//					float3 _dv = {0.f, 0.f, 0.f};
//					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//					{
//						uint j = dNeighbors[dOffsets[i] + _n_idx];
//						float _k_j = density_adv[j] * factor[j];
//						float _k_sum = _k_i + _k_j;
//						if (_k_sum > 1e-5f || _k_sum < -1e-5f)
//						{
//							float3 _grad_p_j = -V[j] * gradW(x[i] - x[j]);
//							_dv -= dt * _k_sum * _grad_p_j;
//						}
//					}
//					v[i] += _dv;
//				});
//		thrust::for_each(
//				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//				[
//						x = Fluid->x.data(),
//						v = Fluid->v.data(),
//						V = Fluid->V.data(),
//						density_adv = Fluid->density_adv.data(),
//						dNeighbors = neighbor_set.d_Neighbors.data(),
//						dCounts = neighbor_set.d_NeighborCounts.data(),
//						dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
//				] __device__(size_t i)
//				{
//					float _d_adv = 0;
//					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//					{
//						uint j = dNeighbors[dOffsets[i] + _n_idx];
//						_d_adv += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
//					}
//					density_adv[i] = max(_d_adv, 0.f);
//				});
//		thrust::transform(Fluid->density_adv.begin(), Fluid->density_adv.end(), Fluid->tmp.begin(), [] __device__(float _) { return _ * REST_DENSITY; });
//		avg_density_error = thrust::reduce(Fluid->tmp.begin(), Fluid->tmp.end(), 0.f, thrust::plus<float>()) / size;
//
//		const float eta = 1.f / dt * 0.1f * 0.01 * REST_DENSITY;
//		if (avg_density_error < eta)
//			break;
//		++iteration_v;
//	}
//	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
//	_) { return _ * dt; });
////	std::cout << "DFSPH - iteration V: " << iteration_v << " Avg density err: " << avg_density_error << std::endl;
//
//
//
//	// ==================== 4. Non-Pressure Force and Predict Velocity ====================
//	thrust::for_each(
//			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//			[
//					x = Fluid->x.data(),
//					v = Fluid->v.data(),
//					a = Fluid->a.data(),
//					m = Fluid->m.data(),
//					rho = Fluid->rho.data(),
//					dNeighbors = neighbor_set.d_Neighbors.data(),
//					dCounts = neighbor_set.d_NeighborCounts.data(),
//					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
//			] __device__(size_t i)
//			{
//				float3 _dv = {0, -9.8f, 0};
//				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//				{
//					uint j = dNeighbors[dOffsets[i] + _n_idx];
//					const float3 _r = x[i] - x[j];
//					const float _r2 = dot(_r, _r);
//					const float _r1 = sqrtf(_r2);
//					const float _diameter = PARTICLE_RADIUS * 2;
//					const float _diameter2 = _diameter * _diameter;
//					if (_r2 > _diameter2)
//						_dv -= SURFACE_TENSION / m[i] * m[j] * _r * W(_r1);
//					else
//						_dv -= SURFACE_TENSION / m[i] * m[j] * _r * W(_diameter);
//
////					float _v_xy = dot(v[i] - v[j], _r);
////					float3 _f_v = 10.f * VISCOSITY * (m[j] / rho[j]) * _v_xy / (_r2 + 0.01f * KERNEL_RADIUS * KERNEL_RADIUS) * gradW(_r);
////					_dv += _f_v;
//				}
//				a[i] = _dv;
//			});
//	thrust::transform(Fluid->a.begin(), Fluid->a.end(), Fluid->v.begin(), Fluid->v.begin(),[dt] __device__(float3
//	a, float3
//	v) { return v + dt * a; });
//
//
//
//	// ==================== 5. Pressure Solve ====================
//	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
//	_) { return _ / (dt * dt); });
//	thrust::for_each(
//			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//			[
//					x = Fluid->x.data(),
//					v = Fluid->v.data(),
//					V = Fluid->V.data(),
//					rho = Fluid->rho.data(),
//					density_adv = Fluid->density_adv.data(),
//					dNeighbors = neighbor_set.d_Neighbors.data(),
//					dCounts = neighbor_set.d_NeighborCounts.data(),
//					dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
//					dt
//			] __device__(size_t i)
//			{
//				float _delta = 0;
//				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//				{
//					uint j = dNeighbors[dOffsets[i] + _n_idx];
//					_delta += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
//				}
//				float _d_adv = rho[i] / REST_DENSITY + dt * _delta;
//				density_adv[i] = max(_d_adv, 1.f);
//			});
//	uint iteration_d = 0;
//	avg_density_error = 0;
//	while (iteration_d < 1 || iteration_d < 100)
//	{
//		thrust::for_each(
//				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//				[
//						x = Fluid->x.data(),
//						v = Fluid->v.data(),
//						V = Fluid->V.data(),
//						factor = Fluid->factor.data(),
//						density_adv = Fluid->density_adv.data(),
//						dNeighbors = neighbor_set.d_Neighbors.data(),
//						dCounts = neighbor_set.d_NeighborCounts.data(),
//						dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
//						dt
//				] __device__(size_t i)
//				{
//					float _k_i = (density_adv[i] - 1.f) * factor[i];
//					float3 _dv = {0.f, 0.f, 0.f};
//					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//					{
//						uint j = dNeighbors[dOffsets[i] + _n_idx];
//						float _k_j = (density_adv[j] - 1.f) * factor[j];
//						float _k_sum = _k_i + _k_j;
//						if (_k_sum > 1e-5f || _k_sum < -1e-5f)
//						{
//							float3 _grad_p_j = -V[j] * gradW(x[i] - x[j]);
//							_dv -= dt * _k_sum * _grad_p_j;
//						}
//					}
//					v[i] += _dv;
//				});
//		thrust::for_each(
//				thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//				[
//						x = Fluid->x.data(),
//						v = Fluid->v.data(),
//						V = Fluid->V.data(),
//						rho = Fluid->rho.data(),
//						density_adv = Fluid->density_adv.data(),
//						dNeighbors = neighbor_set.d_Neighbors.data(),
//						dCounts = neighbor_set.d_NeighborCounts.data(),
//						dOffsets = neighbor_set.d_NeighborWriteOffsets.data(),
//						dt
//				] __device__(size_t i)
//				{
//					float _delta = 0;
//					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
//					{
//						uint j = dNeighbors[dOffsets[i] + _n_idx];
//						_delta += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
//					}
//					float _d_adv = rho[i] / REST_DENSITY + dt * _delta;
//					density_adv[i] = max(_d_adv, 1.f);
//				});
////		avg_density_error = 0;
////		for (int iter = 0; iter < size; ++iter)
////			avg_density_error += REST_DENSITY * (Fluid->density_adv[iter] - 1.f);
////		avg_density_error /= size;
//
//		thrust::transform(Fluid->density_adv.begin(), Fluid->density_adv.end(), Fluid->tmp.begin(), [] __device__(float _) { return (_ - 1.f) * REST_DENSITY; });
//		avg_density_error = thrust::reduce(Fluid->tmp.begin(), Fluid->tmp.end(), 0.f, thrust::plus<float>()) / size;
//
//		const float eta = 0.05f * 0.01f * REST_DENSITY;
//		if (avg_density_error < eta)
//			break;
//		++iteration_d;
//	}
//	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
//	_) { return _ * (dt * dt); });
////	std::cout << "DFSPH - iteration: " << iteration_d << " Avg density err: " << avg_density_error << std::endl;
//
//
//
//	// ==================== 6. Advection ====================
//	thrust::transform(Fluid->v.begin(), Fluid->v.end(), Fluid->x.begin(), Fluid->x.begin(),[dt] __device__(float3
//	v, float3
//	x) { return x + dt * v; });
//
//
//
//	// ==================== 6. Boundary ====================
//	thrust::for_each(
//			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
//			[
//					x = Fluid->x.data(),
//					v = Fluid->v.data()
//			] __device__(size_t i)
//			{
//				float3 collision_normal{0, 0, 0};
//				if (x[i].x > MAX_BOUND.x)
//				{
//					x[i].x = MAX_BOUND.x;
//					collision_normal.x += 1;
//				}
//				if (x[i].x < -MAX_BOUND.x)
//				{
//					x[i].x = -MAX_BOUND.x;
//					collision_normal.x -= 1;
//				}
//				if (!TOP_OPEN)
//				{
//					if (x[i].y > MAX_BOUND.y)
//					{
//						x[i].y = MAX_BOUND.y;
//						collision_normal.y += 1;
//					}
//				}
//				if (x[i].y < -MAX_BOUND.y)
//				{
//					x[i].y = -MAX_BOUND.y;
//					collision_normal.y -= 1;
//				}
//				if (x[i].z > MAX_BOUND.z)
//				{
//					x[i].z = MAX_BOUND.z;
//					collision_normal.z += 1;
//				}
//				if (x[i].z < -MAX_BOUND.z)
//				{
//					x[i].z = -MAX_BOUND.z;
//					collision_normal.z -= 1;
//				}
//				if (length(collision_normal) > std::numeric_limits<float>::epsilon())
//				{
//					collision_normal = normalize(collision_normal);
//					v[i] -= (1. + 0.5f) * dot(v[i], collision_normal) * collision_normal;
//					v[i] *= 0.9f;
//				}
//			});
//}

#ifdef TEST_DFSPH
#include <vector>
#include <numeric>
int main()
{
	using Real = float;
	using Real3 = float3;
	std::vector<Real3> positions;
	std::size_t const N = 120;
	Real const r_omega = static_cast<Real>(0.15);
	Real const r_omega2 = r_omega * r_omega;
	Real const radius = static_cast<Real>(2.0) * (static_cast<Real>(2.0) * r_omega / static_cast<Real>(N - 1));

//Generate test data
	Real min_x = std::numeric_limits<Real>::max();
	Real max_x = std::numeric_limits<Real>::min();
	positions.reserve(N * N * N);
	for (uint i = 0; i < N; ++i)
	{
		for (uint j = 0; j < N; ++j)
		{
			for (uint k = 0; k < N; ++k)
			{
				std::array<Real, 3> x = {{
												 r_omega * static_cast<Real>(2.0 * static_cast<double>(i) / static_cast<double>(N - 1) - 1.0),
												 r_omega * static_cast<Real>(2.0 * static_cast<double>(j) / static_cast<double>(N - 1) - 1.0),
												 r_omega * static_cast<Real>(2.0 * static_cast<double>(k) / static_cast<double>(N - 1) - 1.0)}};

				Real l2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
				if (l2 < r_omega2)
				{
					x[0] += static_cast<Real>(0.35);
					x[1] += static_cast<Real>(0.35);
					x[2] += static_cast<Real>(0.35);
					positions.push_back(make_float3(x[0], x[1], x[2]));
					if (min_x > x[0])
					{
						min_x = x[0];
					}
					if (max_x < x[0])
					{
						max_x = x[0];
					}
				}
			}
		}
	}
	printf("Number of particles: %d \n", static_cast<int>(positions.size()));

	//Create neighborhood search instance
	cuNSearch::NeighborhoodSearch nsearch(radius);

	//Add point set from the test data
	auto pointSetIndex = nsearch.add_point_set(&positions.front().x, positions.size(), true, true);

	for (size_t i = 0; i < 5; i++)
	{
		if (i != 0)
		{
			nsearch.z_sort();
			nsearch.point_set(pointSetIndex).sort_field((Real3 *) nsearch.point_set(pointSetIndex).GetPoints());
		}

		nsearch.find_neighbors();
	}

	//Neighborhood search result test
	auto &pointSet = nsearch.point_set(0);
	auto points = pointSet.GetPoints();

	std::cout << "Validate results" << std::endl;
	for (uint i = 0; i < pointSet.n_points(); i++)
	{
		Real3 point = ((Real3 *) points)[i];
		auto count = pointSet.n_neighbors(0, i);
		for (uint j = 0; j < count; j++)
		{
			auto neighbor = pointSet.neighbor(0, i, j);
			Real3 diff = {point.x - ((Real3 *) points)[neighbor].x, point.y - ((Real3 *) points)[neighbor].y, point.z - ((Real3 *) points)[neighbor].z};
			float squaredLength = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
			float distance = sqrt(squaredLength);

			if (distance > radius)
			{
				throw std::runtime_error("Not a neighbor");
			}
		}
	}
	return 0;
}
#endif
