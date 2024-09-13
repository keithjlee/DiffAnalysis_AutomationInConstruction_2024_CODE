# make constraint function and closure
begin
    function constraint_function(x, p, dmax, fmax)
        res = solve_truss(x, p)

        vertical_displacements = res.U[2:3:end]
        axial_stresses = axial_stress(res, p)

        return [
            (-vertical_displacements .- dmax);
            (axial_stresses .- fmax);
        ]

        # return [
        #     -minimum(vertical_displacements) - dmax,
        #     maximum(axial_stresses) - fmax
        # ]

    end

    CSTR = x -> constraint_function(x, params, dmax, fy)
end

# test constraint function and jacobian
@time c0, dc0 = withjacobian(CSTR, x_init)

# make constraint function and closure
begin
    function constraint_function2(x, p, dmax)
        res = solve_truss(x, p)

        vertical_displacements = res.U[2:3:end]

        return -vertical_displacements .- dmax

        # return [
        #     -minimum(vertical_displacements) - dmax,
        #     maximum(axial_stresses) - fmax
        # ]

    end

    CSTR2 = x -> constraint_function2(x, params, dmax)
end

# test constraint function and jacobian
@time c1, dc1 = withjacobian(CSTR2, x_init)