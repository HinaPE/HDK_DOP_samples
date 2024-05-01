#include "FastMassSpring.h"

float HinaPE::SIMD::SpringConstraint::EvaluateEnergy(const UT_Vector &X) const { return 0; }
void HinaPE::SIMD::SpringConstraint::EvaluateGradient(const UT_Vector &X, UT_Vector &gradient) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateHessian(const UT_Vector &X, std::vector<UT_SparseMatrixCSRF::Triplet> &hessian_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateWeightedLaplacian(std::vector<UT_SparseMatrixCSRF::Triplet> &laplacian_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateWeightedDiagonal(std::vector<UT_SparseMatrixCSRF::Triplet> &diagonal_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateDVector(size_t index, const UT_Vector &X, UT_Vector &d) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateJMatrix(size_t index, std::vector<UT_SparseMatrixCSRF::Triplet> &j_triplets) const {}

HinaPE::SIMD::FastMassSpring::FastMassSpring() : Param()
{
	Cloth = std::make_shared<ClothSIMD>();
}
void HinaPE::SIMD::FastMassSpring::resize(size_t n)
{
	if (size == n)
		return;
	Cloth->x.resize(n, {0, 0, 0});
	Cloth->v.resize(n, {0, 0, 0});
	Cloth->a.resize(n, {0, 0, 0});
	Cloth->m.resize(n, Param.mass);
	Cloth->inv_m.resize(n, 1.f / Param.mass);
	size = n;
}
void HinaPE::SIMD::FastMassSpring::solve_explicit_euler(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_explicit_symplectic(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_implicit_euler(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_gradient_descent(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_newton_descent(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_newton_descent_pcg(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_local_global(float dt) {}

#ifdef TEST_FastMassSpring
#include <iostream>
int main()
{
	std::cout << "Hello, World!" << std::endl;
	return 0;
}
#endif
