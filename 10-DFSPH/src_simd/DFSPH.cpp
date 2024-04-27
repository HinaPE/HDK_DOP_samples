#include "DFSPH.h"
#include "vectorclass.h"

#include <algorithm>

HinaPE::SIMD::DFSPH::DFSPH(float _kernel_radius)
{

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
	Fluid->tmp.resize(n, 0);
	size = n;
}
void HinaPE::SIMD::DFSPH::solve(float dt)
{
	if (size < 16)
		return;

	std::transform(Fluid->a.begin(), Fluid->a.end(), Fluid->a.begin(), [](std::array<float, 3> a) { return std::array<float, 3>{0, -9.8, 0}; });

	size_t left = size % 16;
	for (size_t i = 0; i < size - left; ++i)
	{
		Vec16f x_x(Fluid->x[i][0], Fluid->x[i + 1][0], Fluid->x[i + 2][0], Fluid->x[i + 3][0],
				   Fluid->x[i + 4][0], Fluid->x[i + 5][0], Fluid->x[i + 6][0], Fluid->x[i + 7][0],
				   Fluid->x[i + 8][0], Fluid->x[i + 9][0], Fluid->x[i + 10][0], Fluid->x[i + 11][0],
				   Fluid->x[i + 12][0], Fluid->x[i + 13][0], Fluid->x[i + 14][0], Fluid->x[i + 15][0]);
		Vec16f x_y(Fluid->x[i][1], Fluid->x[i + 1][1], Fluid->x[i + 2][1], Fluid->x[i + 3][1],
				   Fluid->x[i + 4][1], Fluid->x[i + 5][1], Fluid->x[i + 6][1], Fluid->x[i + 7][1],
				   Fluid->x[i + 8][1], Fluid->x[i + 9][1], Fluid->x[i + 10][1], Fluid->x[i + 11][1],
				   Fluid->x[i + 12][1], Fluid->x[i + 13][1], Fluid->x[i + 14][1], Fluid->x[i + 15][1]);
		Vec16f x_z(Fluid->x[i][2], Fluid->x[i + 1][2], Fluid->x[i + 2][2], Fluid->x[i + 3][2],
				   Fluid->x[i + 4][2], Fluid->x[i + 5][2], Fluid->x[i + 6][2], Fluid->x[i + 7][2],
				   Fluid->x[i + 8][2], Fluid->x[i + 9][2], Fluid->x[i + 10][2], Fluid->x[i + 11][2],
				   Fluid->x[i + 12][2], Fluid->x[i + 13][2], Fluid->x[i + 14][2], Fluid->x[i + 15][2]);

		Vec16f v_x(Fluid->v[i][0], Fluid->v[i + 1][0], Fluid->v[i + 2][0], Fluid->v[i + 3][0],
				   Fluid->v[i + 4][0], Fluid->v[i + 5][0], Fluid->v[i + 6][0], Fluid->v[i + 7][0],
				   Fluid->v[i + 8][0], Fluid->v[i + 9][0], Fluid->v[i + 10][0], Fluid->v[i + 11][0],
				   Fluid->v[i + 12][0], Fluid->v[i + 13][0], Fluid->v[i + 14][0], Fluid->v[i + 15][0]);
		Vec16f v_y(Fluid->v[i][1], Fluid->v[i + 1][1], Fluid->v[i + 2][1], Fluid->v[i + 3][1],
				   Fluid->v[i + 4][1], Fluid->v[i + 5][1], Fluid->v[i + 6][1], Fluid->v[i + 7][1],
				   Fluid->v[i + 8][1], Fluid->v[i + 9][1], Fluid->v[i + 10][1], Fluid->v[i + 11][1],
				   Fluid->v[i + 12][1], Fluid->v[i + 13][1], Fluid->v[i + 14][1], Fluid->v[i + 15][1]);
		Vec16f v_z(Fluid->v[i][2], Fluid->v[i + 1][2], Fluid->v[i + 2][2], Fluid->v[i + 3][2],
				   Fluid->v[i + 4][2], Fluid->v[i + 5][2], Fluid->v[i + 6][2], Fluid->v[i + 7][2],
				   Fluid->v[i + 8][2], Fluid->v[i + 9][2], Fluid->v[i + 10][2], Fluid->v[i + 11][2],
				   Fluid->v[i + 12][2], Fluid->v[i + 13][2], Fluid->v[i + 14][2], Fluid->v[i + 15][2]);

		Vec16f a_x(Fluid->a[i][0], Fluid->a[i + 1][0], Fluid->a[i + 2][0], Fluid->a[i + 3][0],
				   Fluid->a[i + 4][0], Fluid->a[i + 5][0], Fluid->a[i + 6][0], Fluid->a[i + 7][0],
				   Fluid->a[i + 8][0], Fluid->a[i + 9][0], Fluid->a[i + 10][0], Fluid->a[i + 11][0],
				   Fluid->a[i + 12][0], Fluid->a[i + 13][0], Fluid->a[i + 14][0], Fluid->a[i + 15][0]);
		Vec16f a_y(Fluid->a[i][1], Fluid->a[i + 1][1], Fluid->a[i + 2][1], Fluid->a[i + 3][1],
				   Fluid->a[i + 4][1], Fluid->a[i + 5][1], Fluid->a[i + 6][1], Fluid->a[i + 7][1],
				   Fluid->a[i + 8][1], Fluid->a[i + 9][1], Fluid->a[i + 10][1], Fluid->a[i + 11][1],
				   Fluid->a[i + 12][1], Fluid->a[i + 13][1], Fluid->a[i + 14][1], Fluid->a[i + 15][1]);
		Vec16f a_z(Fluid->a[i][2], Fluid->a[i + 1][2], Fluid->a[i + 2][2], Fluid->a[i + 3][2],
				   Fluid->a[i + 4][2], Fluid->a[i + 5][2], Fluid->a[i + 6][2], Fluid->a[i + 7][2],
				   Fluid->a[i + 8][2], Fluid->a[i + 9][2], Fluid->a[i + 10][2], Fluid->a[i + 11][2],
				   Fluid->a[i + 12][2], Fluid->a[i + 13][2], Fluid->a[i + 14][2], Fluid->a[i + 15][2]);

		v_x += dt * a_x;
		v_y += dt * a_y;
		v_z += dt * a_z;

		x_x += dt * v_x;
		x_y += dt * v_y;
		x_z += dt * v_z;
	}

	for (size_t i = size - left; i < size; ++i)
	{
		Fluid->v[i][0] += dt * Fluid->a[i][0];
		Fluid->v[i][1] += dt * Fluid->a[i][1];
		Fluid->v[i][2] += dt * Fluid->a[i][2];

		Fluid->x[i][0] += dt * Fluid->v[i][0];
		Fluid->x[i][1] += dt * Fluid->v[i][1];
		Fluid->x[i][2] += dt * Fluid->v[i][2];
	}
}

#ifdef TEST_DFSPH
#include "vectorfp16e.h"
#include <iostream>
int main()
{
	{
		__m128 a = {1.f, 2.f, 3.f, 4.f};
		__m128 b = {5.f, 6.f, 7.f, 8.f};
		__m128 c = _mm_add_ps(a, b);
		__m128 d = _mm_mul_ps(a, b);

		float r_c = _mm_cvtss_f32(c);
		float r_d = _mm_cvtss_f32(d);
		std::cout << r_c << " " << r_d << std::endl;
	}

	{
		__m128i a = {1, 2, 3, 4};
		__m128i b = {5, 6, 7, 8};
		__m128i c = _mm_add_epi32(a, b);
		__m128i d = _mm_mullo_epi32(a, b);
		int r_c = _mm_cvtsi128_si32(c);
		int r_d = _mm_cvtsi128_si32(d);
		std::cout << r_c << " " << r_d << std::endl;
	}

	{
		float a[4] = {1.f, 2.f, 3.f, 4.f};
		__m128 b = _mm_loadu_ps(a);
		__m128 c = _mm_add_ps(b, b);
		float res[4]{};
		_mm_storeu_ps(res, c);
		std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << std::endl;
	}

	{
		Vec16f a(1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8);
		Vec16f b(1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8);
		Vec16f c = a + b;
		std::cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << " " << c[4] << " " << c[5] << " " << c[6] << " " << c[7] << " " << c[8] << " " << c[9] << " " << c[10] << " " << c[11] << " " << c[12] << " " << c[13] << " " << c[14] << " " << c[15] << std::endl;
	}

	return 0;
}
#endif
