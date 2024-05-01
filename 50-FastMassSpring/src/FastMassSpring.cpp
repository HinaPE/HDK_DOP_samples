#include "FastMassSpring.h"

float HinaPE::SIMD::SpringConstraint::EvaluateEnergy(const UT_Vector &X) const { return 0; }
void HinaPE::SIMD::SpringConstraint::EvaluateGradient(const UT_Vector &X, UT_Vector &gradient) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateHessian(const UT_Vector &X, UT_Array<UT_SparseMatrixCSRF::Triplet> &hessian_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateWeightedLaplacian(UT_Array<UT_SparseMatrixCSRF::Triplet> &laplacian_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateWeightedDiagonal(UT_Array<UT_SparseMatrixCSRF::Triplet> &diagonal_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateDVector(size_t index, const UT_Vector &X, UT_Vector &d) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateJMatrix(size_t index, UT_Array<UT_SparseMatrixCSRF::Triplet> &j_triplets) const {}

HinaPE::SIMD::FastMassSpring::FastMassSpring() : Param() { Cloth = std::make_shared<ClothSIMD>(); }
void HinaPE::SIMD::FastMassSpring::build(size_t n, const std::set<std::pair<size_t, size_t>> &springs)
{
	Cloth->x.resize(n, {0, 0, 0});
	Cloth->v.resize(n, {0, 0, 0});
	Cloth->a.resize(n, {0, 0, 0});
	Cloth->f.resize(n, {0, 0, 0});
	Cloth->m.resize(n, Param.mass);
	Cloth->inv_m.resize(n, 1.f / Param.mass);
	Cloth->inertia.resize(n, {0, 0, 0});
	size = n;

	for (const auto &spring: springs)
	{
		SpringConstraint constraint;
		constraint.i = spring.first;
		constraint.j = spring.second;
		constraint.stiffness = Param.stiffness;
		Springs.push_back(constraint);
	}
}
void HinaPE::SIMD::FastMassSpring::solve_explicit_euler(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_explicit_symplectic(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_implicit_euler(float dt)
{
	// compute inertia
	std::transform(Cloth->x.begin(), Cloth->x.end(), Cloth->v.begin(), Cloth->inertia.begin(), [dt](const std::array<float, 3> &x, const std::array<float, 3> &v)
	{
		return std::array<float, 3>{x[0] + dt * v[0], x[1] + dt * v[1], x[2] + dt * v[2]};
	});

	// external force
	std::fill(Cloth->a.begin(), Cloth->a.end(), std::array<float, 3>{0, 0, 0});
	std::transform(Cloth->a.begin(), Cloth->a.end(), Cloth->a.begin(), [this](const std::array<float, 3> &a)
	{
		return std::array<float, 3>{a[0] + Param.gravity[0], a[1] + Param.gravity[1], a[2] + Param.gravity[2]};
	});
	std::transform(Cloth->a.begin(), Cloth->a.end(), Cloth->f.begin(), [this](const std::array<float, 3> &a)
	{
		return std::array<float, 3>{a[0] * Param.mass, a[1] * Param.mass, a[2] * Param.mass};
	});

	// hessian (q_n)
	UT_SparseMatrixCSRF hessian(size * 3, size * 3);
	UT_Array<UT_SparseMatrixCSRF::Triplet> hessian_triplets;
	hessian.setValues(hessian_triplets);
}
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
