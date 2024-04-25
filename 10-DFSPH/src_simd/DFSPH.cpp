#include "DFSPH.h"

HinaPE::SIMD::DFSPH::DFSPH(float _kernel_radius)
{

}
void HinaPE::SIMD::DFSPH::resize(size_t n)
{

}
void HinaPE::SIMD::DFSPH::solve(float dt)
{

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
