import taichi as ti
import taichi.math as tm


@ti.data_oriented
class DFSPHSolver:
    def __init__(self, dt, kernel_radius, particle_radius, rest_density, size, neighbors_max_size):
        # params
        self.dt = dt
        self.kr = kernel_radius
        self.rd = rest_density
        self.sz = size
        self.nsz = neighbors_max_size

        # fields ( init outside )
        self.x = ti.Vector.field(3, dtype=ti.f32, shape=size)
        self.v = ti.Vector.field(3, dtype=ti.f32, shape=size)
        self.a = ti.Vector.field(3, dtype=ti.f32, shape=size)
        self.m = ti.field(dtype=ti.f32, shape=size)
        self.V = ti.field(dtype=ti.f32, shape=size)
        self.rho = ti.field(dtype=ti.f32, shape=size)
        self.factor = ti.field(dtype=ti.f32, shape=size)
        self.neighbors = ti.field(dtype=ti.i32, shape=(size, neighbors_max_size))

        # tmp fields
        self.density_adv = ti.field(dtype=ti.f32, shape=size)

        # fields init
        _V = 0.8 * (2 * particle_radius) ** 3
        _m = self.rd * _V
        self.V.fill(_V)
        self.m.fill(_m)

    @ti.func
    def cubic(self, r_norm):
        res = ti.cast(0.0, ti.f32)
        h = self.kr
        k = 8 / tm.pi
        k /= h ** 3
        q = r_norm / h
        if q <= 1.0:
            if q <= 0.5:
                q2 = q * q
                q3 = q2 * q
                res = k * (6.0 * q3 - 6.0 * q2 + 1)
            else:
                res = k * 2 * ti.pow(1 - q, 3.0)
        return res

    @ti.func
    def d_cubic(self, r):
        h = self.kr
        k = 8 / tm.pi
        k = 6. * k / h ** 3
        r_norm = r.norm()
        q = r_norm / h
        res = ti.Vector([0., 0., 0.])
        if r_norm > 1e-5 and q <= 1.0:
            grad_q = r / (r_norm * h)
            if q <= 0.5:
                res = k * q * (3.0 * q - 2.0) * grad_q
            else:
                factor = 1.0 - q
                res = k * (-factor * factor) * grad_q
        return res

    @ti.kernel
    def compute_density(self):
        for i in ti.grouped(self.x):
            self.rho[i] = self.V[i] * self.cubic(0.0)
            for neighbor in range(self.nsz):
                j = self.neighbors[i, neighbor]
                if j == -1:
                    break
                self.rho[i] += self.V[j] * self.cubic((self.x[i] - self.x[j]).norm())
            self.rho[i] *= self.rd

    @ti.kernel
    def compute_factor(self):
        for i in ti.grouped(self.x):
            sum_grad_p_k = 0.0
            grad_p_i = ti.Vector([0., 0., 0.])
            for neighbor in range(self.nsz):
                j = self.neighbors[i, neighbor]
                if j == -1:
                    break
                grad_p_j = -self.V[j] * self.d_cubic(self.x[i] - self.x[j])
                sum_grad_p_k += grad_p_j.norm_sqr()
                grad_p_i -= grad_p_j
            sum_grad_p_k += grad_p_i.norm_sqr()

            factor = 0.0
            if sum_grad_p_k > 1e-6:
                factor = -1.0 / sum_grad_p_k
            self.factor[i] = factor

    @ti.kernel
    def non_pressure_force(self):
        for i in ti.grouped(self.a):
            self.a[i] = ti.Vector([0, -9.8, 0])

    @ti.kernel
    def predict_velocity(self):
        for i in ti.grouped(self.v):
            self.v[i] += self.dt * self.a[i]

    @ti.kernel
    def advect(self):
        for i in ti.grouped(self.x):
            self.x[i] += self.dt * self.v[i]

    @ti.kernel
    def divergence_solve_kernel(self):
        pass

    def divergence_solve(self):
        pass

    @ti.kernel
    def compute_density_adv(self):
        for i in ti.grouped(self.x):
            delta = 0.
            for neighbor in range(self.nsz):
                j = self.neighbors[i, neighbor]
                if j == -1:
                    break
                delta += self.V[j] * (self.v[i] - self.v[j]).dot(self.d_cubic(self.x[i] - self.x[j]))
            density_adv = self.rho[i] / self.rd + self.dt * delta
            self.density_adv[i] = ti.max(density_adv, 1.0)

    @ti.kernel
    def compute_density_err(self, offset: ti.f32) -> float:
        density_error = 0.0
        for i in ti.grouped(self.x):
            density_error += self.rd * self.density_adv[i] - offset
        return density_error

    @ti.kernel
    def pressure_solve_kernel(self):
        for i in ti.grouped(self.x):
            b_i = self.density_adv[i] - 1.0
            k_i = b_i * self.factor[i] / (self.dt * self.dt)
            for neighbor in range(self.nsz):
                j = self.neighbors[i, neighbor]
                if j == -1:
                    break
                b_j = self.density_adv[j] - 1.0
                k_j = b_j * self.factor[j] / (self.dt * self.dt)
                k_sum = k_i + k_j
                if ti.abs(k_sum) > 1e-5:
                    grad_p_j = -self.V[j] * self.d_cubic(self.x[i] - self.x[j])
                    self.v[i] -= self.dt * k_sum * grad_p_j

    def pressure_solve(self):
        self.compute_density_adv()
        iter = 0
        avg_density_err = 0.0
        while iter < 1 or iter < 100:
            self.pressure_solve_kernel()
            self.compute_density_adv()
            avg_density_err = self.compute_density_err(self.rd) / self.sz
            eta = 0.05 * 0.01 * self.rd
            if ti.abs(avg_density_err) < eta:
                break
            iter += 1
        print(f"DFSPH - iterations: {iter} Avg density Err: {avg_density_err:.4f}")

    def solve(self):
        self.compute_density()
        self.compute_factor()
        self.divergence_solve()
        self.non_pressure_force()
        self.predict_velocity()
        self.pressure_solve()
        self.advect()


if __name__ == '__main__':
    ti.init(arch=ti.gpu)
    solver = DFSPHSolver(0.02, 0.04, 0.01, 1000, 1024, 64)
    solver.solve()
