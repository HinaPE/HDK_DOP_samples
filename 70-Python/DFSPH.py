import taichi as ti
import taichi.math as tm


@ti.data_oriented
class DFSPHSolver:
    def __init__(self, kernel_radius, size, dt):
        # params
        self.dt = dt
        self.kr = kernel_radius
        self.sz = size

        # fields
        self.x = ti.Vector.field(3, dtype=ti.f32, shape=size)
        self.v = ti.Vector.field(3, dtype=ti.f32, shape=size)
        self.a = ti.Vector.field(3, dtype=ti.f32, shape=size)
        self.m = ti.field(dtype=ti.f32, shape=size)
        self.V = ti.field(dtype=ti.f32, shape=size)
        self.rho = ti.field(dtype=ti.f32, shape=size)

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
        # derivative of cubic spline smoothing kernel
        k = 8 / tm.pi
        k = 6. * k / h ** 3
        r_norm = r.norm()
        q = r_norm / h
        res = ti.Vector([0.0 for _ in range(3)])
        if r_norm > 1e-5 and q <= 1.0:
            grad_q = r / (r_norm * h)
            if q <= 0.5:
                res = k * q * (3.0 * q - 2.0) * grad_q
            else:
                factor = 1.0 - q
                res = k * (-factor * factor) * grad_q
        return res

    @ti.kernel
    def mc(self, r_norm: ti.f32) -> float:
        return self.cubic(r_norm)

    @ti.func
    def compute_density_task(self, p_i, p_j, ret: ti.template()):
        x_i = self.x[p_i]
        x_j = self.x[p_j]
        ret += self.V[p_j] * self.cubic((x_i - x_j).norm())

    @ti.kernel
    def compute_density(self):
        for i in ti.grouped(self.x):
            self.rho[i] = 0.0


if __name__ == '__main__':
    ti.init(arch=ti.gpu)
    solver = DFSPHSolver(1.0, 100, 0.01)
    r = solver.mc(0.5)
    print(r)
