#ifndef FAST_MASS_SPRING_H
#define FAST_MASS_SPRING_H

#include <vector>
#include <array>
#include <memory>
#include <set>

#include <UT/UT_SparseMatrix.h>
#include <GU/GU_TriangleMesh.h>

namespace HinaPE::SIMD
{
using ScalarArray = std::vector<float>;
using Vector3Array = std::vector<std::array<float, 3>>;

struct ClothSIMD
{
	Vector3Array x;
	Vector3Array v;
	Vector3Array a;
	Vector3Array f;
	ScalarArray m;
	ScalarArray inv_m;
};

struct ClothParam
{
	float mass = 1.f;
	float stiffness = 1.f;
	std::array<float, 3> gravity;
};

struct SpringConstraint
{
	float EvaluateEnergy(const UT_VectorF &X) const;
	void EvaluateGradient(const UT_VectorF &X, UT_VectorF &gradient) const;
	void EvaluateHessian(const UT_VectorF &X, UT_Array<UT_SparseMatrixCSRF::Triplet> &hessian_triplets) const;
	void EvaluateWeightedLaplacian(UT_Array<UT_SparseMatrixCSRF::Triplet> &laplacian_triplets) const;
	void EvaluateWeightedDiagonal(UT_Array<UT_SparseMatrixCSRF::Triplet> &diagonal_triplets) const;

	void EvaluateDVector(size_t index, const UT_VectorF &X, UT_VectorF &d) const;
	void EvaluateJMatrix(size_t index, UT_Array<UT_SparseMatrixCSRF::Triplet> &j_triplets) const;

	size_t i, j;
	float rest_length;
	float stiffness;
};

struct FastMassSpring
{
	FastMassSpring();
	void build(size_t n, const std::set<std::pair<size_t, size_t>> &springs);
	void solve_explicit_euler(float dt);
	void solve_explicit_symplectic(float dt);
	void solve_implicit_euler(float dt);
	void solve_gradient_descent(float dt);
	void solve_newton_descent(float dt);
	void solve_newton_descent_pcg(float dt);
	void solve_local_global(float dt);

	std::shared_ptr<ClothSIMD> Cloth;
	std::vector<SpringConstraint> Springs;
	ClothParam Param;
private:
	size_t size;
};
}

#endif //FAST_MASS_SPRING_H
