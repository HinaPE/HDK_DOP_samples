# fractal.py
import warp as wp
import numpy as np

wp.init()

num_points = 1024


@wp.kernel
def length(points: wp.array(dtype=wp.vec3),
           lengths: wp.array(dtype=float)):
    # thread index
    tid = wp.tid()
    # compute distance of each point from origin
    lengths[tid] = wp.length(points[tid])


def compute_lengths():
    # allocate an array of 3d points
    points = wp.array(np.random.rand(num_points, 3), dtype=wp.vec3)
    lengths = wp.zeros(num_points, dtype=float)

    # launch kernel
    wp.launch(kernel=length,
              dim=len(points),
              inputs=[points, lengths])

    return lengths


if __name__ == '__main__':
    lengths = compute_lengths()
    print(lengths)
