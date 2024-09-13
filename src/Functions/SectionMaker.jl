#tube sizing
a_tube(d, alpha) = .25pi * d^2 * (1 - alpha^2)
I_tube(d, alpha) = pi / 64 * d^4 * (1 - alpha^4)
S_tube(d, alpha) = 2 * I_tube(d, alpha) / d
J_tube(d, alpha) = 2 * I_tube(d, alpha)

a_box(b, d) = b * d
Ix_box(b, d) = b * d^3 / 12
Iy_box(b, d) = d * b^3 / 12
Sx_box(b, d) = 2 * Ix_box(b, d) / d
Sy_box(b, d) = 2 * Iy_box(b, d) / b

import Interpolations
ab_ratios = [1., 1.5, 2., 2.5, 3., 4., 5., 6., 10.]
β = [.141, .196, .229, .249, .263, .281, .291, .299, .312]
J_interpolator = Interpolations.linear_interpolation(ab_ratios, β)
function J_box(width, depth)
    ab = max(depth/width, width/depth)

    if ab > 10
        return 1 / 3
    else
        return J_interpolator(ab) * max(depth, width) * min(depth, width)^3
    end

end
