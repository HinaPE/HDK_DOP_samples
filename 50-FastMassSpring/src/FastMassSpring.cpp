#include "FastMassSpring.h"

float HinaPE::SIMD::SpringConstraint::EvaluateEnergy(const UT_VectorF &X) const
{
	UT_Vector3F x_i = {X(i * 3), X(i * 3 + 1), X(i * 3 + 2)};
	UT_Vector3F x_j = {X(j * 3), X(j * 3 + 1), X(j * 3 + 2)};
	UT_Vector3F x_ij = x_i - x_j;
	float d = x_ij.length() - rest_length;
	return 0.5f * stiffness * d * d;
}
void HinaPE::SIMD::SpringConstraint::EvaluateGradient(const UT_VectorF &X, UT_VectorF &gradient) const
{
	UT_Vector3F x_i = {X(i * 3), X(i * 3 + 1), X(i * 3 + 2)};
	UT_Vector3F x_j = {X(j * 3), X(j * 3 + 1), X(j * 3 + 2)};
	UT_Vector3F x_ij = x_i - x_j;
	UT_Vector3F x_normalized = x_ij;
	x_normalized.normalize();
	UT_Vector3F g_ij = stiffness * (x_ij.length() - rest_length) * x_normalized;
	gradient(i * 3) += g_ij.x();
	gradient(i * 3 + 1) += g_ij.y();
	gradient(i * 3 + 2) += g_ij.z();
	gradient(j * 3) -= g_ij.x();
	gradient(j * 3 + 1) -= g_ij.y();
	gradient(j * 3 + 2) -= g_ij.z();
}
void HinaPE::SIMD::SpringConstraint::EvaluateHessian(const UT_VectorF &X, UT_Array<UT_SparseMatrixCSRF::Triplet> &hessian_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateWeightedLaplacian(UT_Array<UT_SparseMatrixCSRF::Triplet> &laplacian_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateWeightedDiagonal(UT_Array<UT_SparseMatrixCSRF::Triplet> &diagonal_triplets) const {}
void HinaPE::SIMD::SpringConstraint::EvaluateDVector(size_t index, const UT_VectorF &X, UT_VectorF &d) const {}
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
	UT_Vector inertia;
	UT_Vector b;

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
