#ifndef FAST_MASS_SPRING_H
#define FAST_MASS_SPRING_H

#include <vector>
#include <array>
#include <memory>

#include <UT/UT_SparseMatrix.h>

namespace HinaPE::SIMD
{
using ScalarArray = std::vector<float>;
using Vector3Array = std::vector<std::array<float, 3>>;

struct ClothSIMD
{
	Vector3Array x;
	Vector3Array v;
	Vector3Array a;
	ScalarArray m;
	ScalarArray inv_m;
};

struct ClothParam
{
	float mass = 1.f;
	std::array<float, 3> gravity;
};

struct SpringConstraint
{
	float EvaluateEnergy(const UT_Vector &X) const;
	void EvaluateGradient(const UT_Vector &X, UT_Vector &gradient) const;
	void EvaluateHessian(const UT_Vector &X, std::vector<UT_SparseMatrixCSRF::Triplet> &hessian_triplets) const;
	void EvaluateWeightedLaplacian(std::vector<UT_SparseMatrixCSRF::Triplet> &laplacian_triplets) const;
	void EvaluateWeightedDiagonal(std::vector<UT_SparseMatrixCSRF::Triplet> &diagonal_triplets) const;

	void EvaluateDVector(size_t index, const UT_Vector &X, UT_Vector &d) const;
	void EvaluateJMatrix(size_t index, std::vector<UT_SparseMatrixCSRF::Triplet> &j_triplets) const;

	float stiffness;
};

struct FastMassSpring
{
	FastMassSpring();
	void resize(size_t n);
	void solve_explicit_euler(float dt);
	void solve_explicit_symplectic(float dt);
	void solve_implicit_euler(float dt);
	void solve_gradient_descent(float dt);
	void solve_newton_descent(float dt);
	void solve_newton_descent_pcg(float dt);
	void solve_local_global(float dt);

	std::shared_ptr<ClothSIMD> Cloth;
	ClothParam Param;
private:
	size_t size;
};
}

#endif //FAST_MASS_SPRING_H
