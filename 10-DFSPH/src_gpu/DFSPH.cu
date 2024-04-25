#include "DFSPH.cuh"
#include <cuda_runtime.h>
#include "common/helper_cuda.h"
#include "neighbor/include/cuNSearch.h"
#include "neighbor/src/cuNSearchDeviceData.h"

__constant__ float KERNEL_RADIUS;
__constant__ float KERNEL_K;
__constant__ float KERNEL_L;
static inline __device__ float dot(const float3 a, const float3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline __device__ float length(const float3 r) { return sqrtf(dot(r, r)); }
static inline __device__ float W(const float3 r)
{
	float res = 0.f;
	const float q = length(r) / KERNEL_RADIUS;
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
			res.x = KERNEL_L * q * (3.f * q - 2.f) * gradq.x;
			res.y = KERNEL_L * q * (3.f * q - 2.f) * gradq.y;
			res.z = KERNEL_L * q * (3.f * q - 2.f) * gradq.z;
		} else
		{
			const float factor = 1.f - q;
			res.x = KERNEL_L * (-factor * factor) * gradq.x;
			res.y = KERNEL_L * (-factor * factor) * gradq.y;
			res.z = KERNEL_L * (-factor * factor) * gradq.z;
		}
	}
	return res;
}
static inline __device__ float W0() { return W(make_float3(0.f, 0.f, 0.f)); }

template<class Func>
__global__ void parallel_for(int n, Func func) { for (int i = blockDim.x * blockIdx.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) func(i); }

HinaPE::CUDA::DFSPH::DFSPH(float _kernel_radius) : kernel_radius(_kernel_radius), size(0)
{
	float k = 8.f / (3.14159265358979323846f * kernel_radius * kernel_radius * kernel_radius);
	float l = 48.f / (3.14159265358979323846f * kernel_radius * kernel_radius * kernel_radius);
	cudaMemcpyToSymbol(KERNEL_RADIUS, &kernel_radius, sizeof(float));
	cudaMemcpyToSymbol(KERNEL_K, &k, sizeof(float));
	cudaMemcpyToSymbol(KERNEL_L, &l, sizeof(float));

	Fluid = std::make_shared<FluidGPU>();
	Searcher = std::make_shared<cuNSearch::NeighborhoodSearch>(kernel_radius);
	fluid_idx = Searcher->add_point_set(&(Fluid->x.data()->x), Fluid->x.size(), true, true, true);
}

void HinaPE::CUDA::DFSPH::resize(size_t n)
{
	Fluid->x.resize(n);
	Fluid->v.resize(n);
	Fluid->a.resize(n);
	Fluid->m.resize(n);
	Fluid->V.resize(n);
	Fluid->rho.resize(n);
	size = n;
}

void HinaPE::CUDA::DFSPH::solve(float dt)
{
	Searcher->update_point_set(fluid_idx);
	Searcher->find_neighbors();
	thrust::for_each(
			thrust::make_counting_iterator((size_t) 0), thrust::make_counting_iterator(size),
			[
					x = Fluid->x.data(),
					v = Fluid->v.data(),
					a = Fluid->a.data(),
					m = Fluid->m.data(),
					V = Fluid->V.data(),
					rho = Fluid->rho.data()
			] __device__(size_t i)
			{
				// compute density
				float rho_i = V[i] * W0();
			});
}

#ifdef TEST_DFSPH
int main()
{
	int n = 65536;
	HinaPE::CUDA::DFSPH df(0.04f);
	df.resize(n);
	df.solve(0.02f);
	return 0;
}
#endif
