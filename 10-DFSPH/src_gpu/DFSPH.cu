#include "DFSPH.cuh"
#include <cuda_runtime.h>
#include "common/helper_cuda.h"
#include "common/helper_math.h"
#include "neighbor/include/cuNSearch.h"
#include "neighbor/src/cuNSearchDeviceData.h"
#include <thrust/execution_policy.h>

__constant__ float PARTICLE_RADIUS;
__constant__ float KERNEL_RADIUS;
__constant__ float KERNEL_K;
__constant__ float KERNEL_L;
__constant__ float REST_DENSITY;
__constant__ float VISCOSITY;
__constant__ float SURFACE_TENSION;
__constant__ float3 MAX_BOUND;
__constant__ bool TOP_OPEN;
static inline __device__ float W(const float r)
{
	float res = 0.f;
	const float q = r / KERNEL_RADIUS;
	if (q <= 1.f)
	{
		if (q <= 0.5f)
		{
			const float q2 = q * q;
			const float q3 = q2 * q;
			res = KERNEL_K * (6.f * q3 - 6.f * q2 + 1.f);
		} else
		{
			res = KERNEL_K * (2.f * powf(1.f - q, 3.f));
		}
	}
	return res;
}
static inline __device__ float W(const float3 r) { return W(length(r)); }
static inline __device__ float3 gradW(const float3 r)
{
	float3 res = {0.f, 0.f, 0.f};
	const float rl = length(r);
	const float q = rl / KERNEL_RADIUS;
	if ((rl > 1.0e-9) && (q <= 1.f))
	{
		float3 gradq = {r.x / rl, r.y / rl, r.z / rl};
		gradq.x /= KERNEL_RADIUS;
		gradq.y /= KERNEL_RADIUS;
		gradq.z /= KERNEL_RADIUS;
		if (q <= 0.5f)
		{
			res = KERNEL_L * q * (3.f * q - 2.f) * gradq;
		} else
		{
			const float factor = 1.f - q;
			res = KERNEL_L * (-factor * factor) * gradq;
		}
	}
	return res;
}
static inline __device__ float W0() { return W(make_float3(0.f, 0.f, 0.f)); }

HinaPE::CUDA::DFSPH::DFSPH(float _kernel_radius) : size(0)
{
	float kr = _kernel_radius;
	float r = 0.01f;
	float k = 8.f / (3.14159265358979323846f * kr * kr * kr);
	float l = 48.f / (3.14159265358979323846f * kr * kr * kr);
	float rd = 1000.f;
	float vis = 0.01f;
	float st = 0.01f;
	float3 max_bound = {2.f, 2.f, 2.f};
	bool top_open = true;
	cudaMalloc((void **) &PARTICLE_RADIUS, sizeof(float));
	cudaMalloc((void **) &KERNEL_RADIUS, sizeof(float));
	cudaMalloc((void **) &KERNEL_K, sizeof(float));
	cudaMalloc((void **) &KERNEL_L, sizeof(float));
	cudaMalloc((void **) &REST_DENSITY, sizeof(float));
	cudaMalloc((void **) &VISCOSITY, sizeof(float));
	cudaMalloc((void **) &SURFACE_TENSION, sizeof(float));
	cudaMalloc((void **) &MAX_BOUND, sizeof(float3));
	cudaMalloc((void **) &TOP_OPEN, sizeof(bool));
	cudaMemcpyToSymbol(PARTICLE_RADIUS, &r, sizeof(float));
	cudaMemcpyToSymbol(KERNEL_RADIUS, &kr, sizeof(float));
	cudaMemcpyToSymbol(KERNEL_K, &k, sizeof(float));
	cudaMemcpyToSymbol(KERNEL_L, &l, sizeof(float));
	cudaMemcpyToSymbol(REST_DENSITY, &rd, sizeof(float));
	cudaMemcpyToSymbol(VISCOSITY, &vis, sizeof(float));
	cudaMemcpyToSymbol(SURFACE_TENSION, &st, sizeof(float));
	cudaMemcpyToSymbol(MAX_BOUND, &max_bound, 3 * sizeof(float));
	cudaMemcpyToSymbol(TOP_OPEN, &top_open, sizeof(bool));

	Fluid = std::make_shared<FluidGPU>();
	Searcher = std::make_shared<cuNSearch::NeighborhoodSearch>(kr);
	fluid_idx = Searcher->add_point_set(&(Fluid->x.data()->x), Fluid->x.size(), true, true, true);
	Searcher->set_active(true);
}

void HinaPE::CUDA::DFSPH::resize(size_t n)
{
	if (size == n)
		return;
	Fluid->x.resize(n, make_float3(0, 0, 0));
	Fluid->v.resize(n, make_float3(0, 0, 0));
	Fluid->a.resize(n, make_float3(0, 0, 0));
	Fluid->m.resize(n, 1000.f * 0.9f * 0.02f * 0.02f * 0.02f);
	Fluid->V.resize(n, 0.9f * 0.02f * 0.02f * 0.02f);
	Fluid->rho.resize(n, 0);
	Fluid->factor.resize(n, 0);
	Fluid->density_adv.resize(n, 0);
	Fluid->nn.resize(n, 0);
	Fluid->tmp.resize(n, 0);
	size = n;

	need_reload = true;
}

void HinaPE::CUDA::DFSPH::solve(float dt)
{

	// ==================== 1. Build Neighbors ====================
	if (need_reload)
	{
		Searcher->resize_point_set(fluid_idx, &(Fluid->x.data()->x), size);
		need_reload = false;
	}
	Searcher->update_point_set(fluid_idx);
	Searcher->find_neighbors();
	cuNSearch::PointSet::NeighborSet &neighbor_set = Searcher->point_set(fluid_idx).get_raw_neighbor_set(fluid_idx);



	// ==================== 2. Compute Density and Factor ====================
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
				float _rho_i = V[i] * W0();
				float _sum_grad_p_k = 0;
				float3 _grad_p_i = {0.f, 0.f, 0.f};
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];
					_rho_i += V[j] * W(x[i] - x[j]);

					float3 _grad_p_j = -V[j] * gradW(x[i] - x[j]);
					_sum_grad_p_k += dot(_grad_p_j, _grad_p_j);
					_grad_p_i -= _grad_p_j;
				}

				_sum_grad_p_k += dot(_grad_p_i, _grad_p_i);
				rho[i] = _rho_i * REST_DENSITY;
				if (_sum_grad_p_k > 1e-6f)
					factor[i] = -1.f / _sum_grad_p_k;
				else
					factor[i] = 0;

				nn[i] = dCounts[i];
			});



	// ==================== 3. Divergence Solve ====================
	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float _) { return _ / dt; });
	thrust::for_each(
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
				float _d_adv = 0;
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];
					_d_adv += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
				}
				density_adv[i] = max(_d_adv, 0.f);
			});
	uint iteration_v = 0;
	float avg_density_error = 0;
	while (iteration_v < 1 || iteration_v < 100)
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
					float _k_i = density_adv[i] * factor[i];
					float3 _dv = {0.f, 0.f, 0.f};
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						float _k_j = density_adv[j] * factor[j];
						float _k_sum = _k_i + _k_j;
						if (_k_sum > 1e-5f || _k_sum < -1e-5f)
						{
							float3 _grad_p_j = -V[j] * gradW(x[i] - x[j]);
							_dv -= dt * _k_sum * _grad_p_j;
						}
					}
					v[i] += _dv;
				});
		thrust::for_each(
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
					float _d_adv = 0;
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						_d_adv += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
					}
					density_adv[i] = max(_d_adv, 0.f);
				});
		thrust::transform(Fluid->density_adv.begin(), Fluid->density_adv.end(), Fluid->tmp.begin(), [] __device__(float _) { return _ * REST_DENSITY; });
		avg_density_error = thrust::reduce(Fluid->tmp.begin(), Fluid->tmp.end(), 0.f, thrust::plus<float>()) / size;

		const float eta = 1.f / dt * 0.1f * 0.01 * REST_DENSITY;
		if (avg_density_error < eta)
			break;
		++iteration_v;
	}
	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
	_) { return _ * dt; });
//	std::cout << "DFSPH - iteration V: " << iteration_v << " Avg density err: " << avg_density_error << std::endl;



	// ==================== 4. Non-Pressure Force and Predict Velocity ====================
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data(),
					a = Fluid->a.data(),
					m = Fluid->m.data(),
					rho = Fluid->rho.data(),
					dNeighbors = neighbor_set.d_Neighbors.data(),
					dCounts = neighbor_set.d_NeighborCounts.data(),
					dOffsets = neighbor_set.d_NeighborWriteOffsets.data()
			] __device__(size_t i)
			{
				float3 _dv = {0, -9.8f, 0};
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];
					const float3 _r = x[i] - x[j];
					const float _r2 = dot(_r, _r);
					const float _r1 = sqrtf(_r2);
					const float _diameter = PARTICLE_RADIUS * 2;
					const float _diameter2 = _diameter * _diameter;
					if (_r2 > _diameter2)
						_dv -= SURFACE_TENSION / m[i] * m[j] * _r * W(_r1);
					else
						_dv -= SURFACE_TENSION / m[i] * m[j] * _r * W(_diameter);

//					float _v_xy = dot(v[i] - v[j], _r);
//					float3 _f_v = 10.f * VISCOSITY * (m[j] / rho[j]) * _v_xy / (_r2 + 0.01f * KERNEL_RADIUS * KERNEL_RADIUS) * gradW(_r);
//					_dv += _f_v;
				}
				a[i] = _dv;
			});
	thrust::transform(Fluid->a.begin(), Fluid->a.end(), Fluid->v.begin(), Fluid->v.begin(),[dt] __device__(float3
	a, float3
	v) { return v + dt * a; });



	// ==================== 5. Pressure Solve ====================
	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
	_) { return _ / (dt * dt); });
	thrust::for_each(
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
				float _delta = 0;
				for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
				{
					uint j = dNeighbors[dOffsets[i] + _n_idx];
					_delta += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
				}
				float _d_adv = rho[i] / REST_DENSITY + dt * _delta;
				density_adv[i] = max(_d_adv, 1.f);
			});
	uint iteration_d = 0;
	avg_density_error = 0;
	while (iteration_d < 1 || iteration_d < 100)
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
					float _k_i = (density_adv[i] - 1.f) * factor[i];
					float3 _dv = {0.f, 0.f, 0.f};
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						float _k_j = (density_adv[j] - 1.f) * factor[j];
						float _k_sum = _k_i + _k_j;
						if (_k_sum > 1e-5f || _k_sum < -1e-5f)
						{
							float3 _grad_p_j = -V[j] * gradW(x[i] - x[j]);
							_dv -= dt * _k_sum * _grad_p_j;
						}
					}
					v[i] += _dv;
				});
		thrust::for_each(
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
					float _delta = 0;
					for (uint _n_idx = 0; _n_idx < dCounts[i]; ++_n_idx)
					{
						uint j = dNeighbors[dOffsets[i] + _n_idx];
						_delta += V[j] * dot(v[i] - v[j], gradW(x[i] - x[j]));
					}
					float _d_adv = rho[i] / REST_DENSITY + dt * _delta;
					density_adv[i] = max(_d_adv, 1.f);
				});
//		avg_density_error = 0;
//		for (int iter = 0; iter < size; ++iter)
//			avg_density_error += REST_DENSITY * (Fluid->density_adv[iter] - 1.f);
//		avg_density_error /= size;

		thrust::transform(Fluid->density_adv.begin(), Fluid->density_adv.end(), Fluid->tmp.begin(), [] __device__(float _) { return (_ - 1.f) * REST_DENSITY; });
		avg_density_error = thrust::reduce(Fluid->tmp.begin(), Fluid->tmp.end(), 0.f, thrust::plus<float>()) / size;

		const float eta = 0.05f * 0.01f * REST_DENSITY;
		if (avg_density_error < eta)
			break;
		++iteration_d;
	}
	thrust::transform(Fluid->factor.begin(), Fluid->factor.end(), Fluid->factor.begin(),[dt] __device__(float
	_) { return _ * (dt * dt); });
//	std::cout << "DFSPH - iteration: " << iteration_d << " Avg density err: " << avg_density_error << std::endl;



	// ==================== 6. Advection ====================
	thrust::transform(Fluid->v.begin(), Fluid->v.end(), Fluid->x.begin(), Fluid->x.begin(),[dt] __device__(float3 v, float3 x) { return x + dt * v; });



	// ==================== 6. Boundary ====================
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data()
			] __device__(size_t i)
			{
				float3 collision_normal{0, 0, 0};
				if (x[i].x > MAX_BOUND.x)
				{
					x[i].x = MAX_BOUND.x;
					collision_normal.x += 1;
				}
				if (x[i].x < -MAX_BOUND.x)
				{
					x[i].x = -MAX_BOUND.x;
					collision_normal.x -= 1;
				}
				if (!TOP_OPEN)
				{
					if (x[i].y > MAX_BOUND.y)
					{
						x[i].y = MAX_BOUND.y;
						collision_normal.y += 1;
					}
				}
				if (x[i].y < -MAX_BOUND.y)
				{
					x[i].y = -MAX_BOUND.y;
					collision_normal.y -= 1;
				}
				if (x[i].z > MAX_BOUND.z)
				{
					x[i].z = MAX_BOUND.z;
					collision_normal.z += 1;
				}
				if (x[i].z < -MAX_BOUND.z)
				{
					x[i].z = -MAX_BOUND.z;
					collision_normal.z -= 1;
				}
				if (length(collision_normal) > std::numeric_limits<float>::epsilon())
				{
					collision_normal = normalize(collision_normal);
					v[i] -= (1. + 0.5f) * dot(v[i], collision_normal) * collision_normal;
					v[i] *= 0.9f;
				}
			});
}

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
