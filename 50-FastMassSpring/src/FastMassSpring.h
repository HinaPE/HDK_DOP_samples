#ifndef FAST_MASS_SPRING_H
#define FAST_MASS_SPRING_H

#include <vector>
#include <array>
#include <memory>
#include <set>

#include <UT/UT_SparseMatrix.h>
#include <GU/GU_TriangleMesh.h>
#include <UT/UT_Vector.h>
#include <UT/UT_Vector3.h>

namespace HinaPE::SIMD
{
struct ClothSIMD
{
	UT_VectorF x;
	UT_VectorF v;
	UT_VectorF f;
	UT_SparseMatrixCSRF m;
	UT_SparseMatrixCSRF inv_m;

	// temp variables
	UT_VectorF y;
	UT_VectorF x_next;
	UT_VectorF v_next;
};

struct ClothParam
{
	float mass = 1.f;
	float stiffness = 80.f;
	std::array<float, 3> gravity;
};

struct AttachmentConstraint
{
	size_t i;
	std::array<float, 3> p;
};

struct SpringConstraint
{
	size_t i, j;
	float rest_length;
	float stiffness;
};

struct FastMassSpring
{
	FastMassSpring();
	void build(const std::vector<std::array<float, 3>> &init_x, const std::set<std::pair<size_t, size_t>> &springs);
	void solve_explicit_euler(float dt);
	void solve_explicit_symplectic(float dt);
	void solve_implicit_euler(float dt);
	void solve_gradient_descent(float dt);
	void solve_newton_descent(float dt);
	void solve_newton_descent_pcg(float dt);
	void solve_local_global(float dt);

	std::shared_ptr<ClothSIMD> Cloth;
	std::vector<AttachmentConstraint> Attachments;
	std::vector<SpringConstraint> Springs;
	ClothParam Param;
private:
	size_t size;
};
}

#endif //FAST_MASS_SPRING_H
