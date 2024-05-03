#include "FastMassSpring.h"

HinaPE::SIMD::FastMassSpring::FastMassSpring() : Param() { Cloth = std::make_shared<ClothSIMD>(); }
void HinaPE::SIMD::FastMassSpring::build(const std::vector<std::array<float, 3>> &init_x, const std::set<std::pair<size_t, size_t>> &springs)
{
	size = init_x.size();
	exint NL = 0;
	exint NH = 3 * size - 1;
	Cloth->x.init(NL, NH);
	Cloth->v.init(NL, NH);
	Cloth->f.init(NL, NH);
	Cloth->m.init(3 * size, 3 * size);
	Cloth->inv_m.init(3 * size, 3 * size);
	Cloth->y.init(NL, NH);
	Cloth->tmp.init(NL, NH);
	Cloth->x_next.init(NL, NH);
	Cloth->v_next.init(NL, NH);

	Cloth->x.zero();
	Cloth->v.zero();
	Cloth->f.zero();
	Cloth->m.zero();
	Cloth->inv_m.zero();
	Cloth->y.zero();
	Cloth->tmp.zero();
	Cloth->x_next.zero();
	Cloth->v_next.zero();
	{
		UT_Array<UT_SparseMatrixCSRF::Triplet> m_triplets;
		UT_Array<UT_SparseMatrixCSRF::Triplet> inv_m_triplets;
		for (size_t i = 0; i < size; ++i)
		{
			m_triplets.append(UT_SparseMatrixCSRF::Triplet(3 * i + 0, 3 * i + 0, Param.mass));
			m_triplets.append(UT_SparseMatrixCSRF::Triplet(3 * i + 1, 3 * i + 1, Param.mass));
			m_triplets.append(UT_SparseMatrixCSRF::Triplet(3 * i + 2, 3 * i + 2, Param.mass));
			inv_m_triplets.append(UT_SparseMatrixCSRF::Triplet(3 * i + 0, 3 * i + 0, 1.f / Param.mass));
			inv_m_triplets.append(UT_SparseMatrixCSRF::Triplet(3 * i + 1, 3 * i + 1, 1.f / Param.mass));
			inv_m_triplets.append(UT_SparseMatrixCSRF::Triplet(3 * i + 2, 3 * i + 2, 1.f / Param.mass));
		}
		Cloth->m.setValues(m_triplets);
		Cloth->inv_m.setValues(inv_m_triplets);
	}

	Attachments.clear();
	Attachments.emplace_back(AttachmentConstraint{0, {init_x[0][0], init_x[0][1], init_x[0][2]}});

	Springs.clear();
	for (const auto &spring: springs)
	{
		SpringConstraint constraint;
		constraint.i = spring.first;
		constraint.j = spring.second;
		constraint.rest_length = std::sqrt(std::pow(init_x[spring.first][0] - init_x[spring.second][0], 2) +
										   std::pow(init_x[spring.first][1] - init_x[spring.second][1], 2) +
										   std::pow(init_x[spring.first][2] - init_x[spring.second][2], 2));
		constraint.stiffness = Param.stiffness;
		Springs.push_back(constraint);
	}
}
void HinaPE::SIMD::FastMassSpring::solve_explicit_euler(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_explicit_symplectic(float dt)
{
	exint NL = 0;
	exint NH = 3 * size - 1;

	for (int i = 0; i < Cloth->f.length(); i += 3)
	{
		Cloth->f(i + 0) = Param.gravity[0];
		Cloth->f(i + 1) = Param.gravity[1];
		Cloth->f(i + 2) = Param.gravity[2];
	}

	for (const auto &attachment: Attachments)
	{
		UT_Vector3F g_i = Param.stiffness * UT_Vector3F(Cloth->x(attachment.i * 3 + 0) - attachment.p[0],
														Cloth->x(attachment.i * 3 + 1) - attachment.p[1],
														Cloth->x(attachment.i * 3 + 2) - attachment.p[2]);
		Cloth->f(attachment.i * 3 + 0) -= g_i.x();
		Cloth->f(attachment.i * 3 + 1) -= g_i.y();
		Cloth->f(attachment.i * 3 + 2) -= g_i.z();
	}
	for (const auto &spring: Springs)
	{
		UT_Vector3F x_ij = UT_Vector3F(Cloth->x(spring.i * 3 + 0) - Cloth->x(spring.j * 3 + 0),
									   Cloth->x(spring.i * 3 + 1) - Cloth->x(spring.j * 3 + 1),
									   Cloth->x(spring.i * 3 + 2) - Cloth->x(spring.j * 3 + 2));
		UT_Vector3F dir_ij = x_ij;
		dir_ij.normalize();
		UT_Vector3F grad_ij = spring.stiffness * (x_ij.length() - spring.rest_length) * dir_ij;
		Cloth->f(spring.i * 3 + 0) -= grad_ij.x();
		Cloth->f(spring.i * 3 + 1) -= grad_ij.y();
		Cloth->f(spring.i * 3 + 2) -= grad_ij.z();
		Cloth->f(spring.j * 3 + 0) += grad_ij.x();
		Cloth->f(spring.j * 3 + 1) += grad_ij.y();
		Cloth->f(spring.j * 3 + 2) += grad_ij.z();
	}

	Cloth->tmp = Cloth->v;
	Cloth->tmp *= dt;
	Cloth->y = Cloth->x;
	Cloth->y += Cloth->tmp;

	Cloth->inv_m.multVec(Cloth->f, Cloth->x_next);
	Cloth->x_next *= dt * dt;
	Cloth->x_next += Cloth->y;

	Cloth->v_next = Cloth->x_next;
	Cloth->v_next -= Cloth->x;
	Cloth->v_next /= dt;

	Cloth->x = Cloth->x_next;
	Cloth->v = Cloth->v_next;
}
void HinaPE::SIMD::FastMassSpring::solve_implicit_euler(float dt)
{
//	// compute inertia
//	UT_VectorF inertia;
//	UT_VectorF b;
//
//	// external force
//	std::fill(Cloth->a.begin(), Cloth->a.end(), std::array<float, 3>{0, 0, 0});
//	std::transform(Cloth->a.begin(), Cloth->a.end(), Cloth->a.begin(), [this](const std::array<float, 3> &a)
//	{
//		return std::array<float, 3>{a[0] + Param.gravity[0], a[1] + Param.gravity[1], a[2] + Param.gravity[2]};
//	});
//	std::transform(Cloth->a.begin(), Cloth->a.end(), Cloth->f.begin(), [this](const std::array<float, 3> &a)
//	{
//		return std::array<float, 3>{a[0] * Param.mass, a[1] * Param.mass, a[2] * Param.mass};
//	});
//
//	// hessian (q_n)
//	UT_SparseMatrixCSRF hessian(size * 3, size * 3);
//	UT_Array<UT_SparseMatrixCSRF::Triplet> hessian_triplets;
//	for (const auto &spring: Springs)
//	{
//		spring.EvaluateHessian(UT_VectorF(), hessian_triplets);
//	}
//
//	hessian.multVec(b, b);
}
void HinaPE::SIMD::FastMassSpring::solve_gradient_descent(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_newton_descent(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_newton_descent_pcg(float dt) {}
void HinaPE::SIMD::FastMassSpring::solve_local_global(float dt) {}

#ifdef TEST_FastMassSpring
#include <iostream>
int main()
{
//	UT_SparseMatrixCSRF A(2, 2);
//	UT_Array<UT_SparseMatrixCSRF::Triplet> A_triplets;
//	A_triplets.append(UT_SparseMatrixCSRF::Triplet(0, 0, 4));
//	A_triplets.append(UT_SparseMatrixCSRF::Triplet(0, 1, 12));
//	A_triplets.append(UT_SparseMatrixCSRF::Triplet(1, 0, 12));
//	A_triplets.append(UT_SparseMatrixCSRF::Triplet(1, 1, 37));
//	A.setValues(A_triplets);
//	UT_VectorF b(0, 1);
//	b(0) = 16;
//	b(1) = 43;
//
//	UT_SparseMatrixCSRF M(2, 2);
//	UT_SparseMatrixCSRF MT(2, 2);
//	UT_SparseMatrixCSRF R(2, 2);
//	M = A;
//	printf("Result: %d\n", M.incompleteCholeskyFactorization());
//
//	M.transpose(MT);
//	MT.multMatrix(M, R);
//	R.printFull(std::cout);
//
//	UT_VectorF x(0, 1);
//	UT_VectorF y(0, 1);
//	MT.solveLowerTriangular(y, b);
//	M.solveUpperTriangular(x, y);
//	std::cout << x << std::endl;


	UT_VectorF F(0, 3 * 100 - 1);
	UT_VectorF _external_force(0, 3 * 100 - 1);
	{
		for (int i = 0; i < _external_force.length(); i += 3)
		{
			_external_force(i + 0) = 0;
			_external_force(i + 1) = -9.8;
			_external_force(i + 2) = 0;
		}
	}
	F = _external_force;
	std::cout << _external_force << std::endl;
	std::cout << F << std::endl;

	return 0;
}
#endif