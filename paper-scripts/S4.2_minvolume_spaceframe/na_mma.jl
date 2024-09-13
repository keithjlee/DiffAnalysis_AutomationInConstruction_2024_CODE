include("init_problem.jl")

#=
OPTIMIZATION: FUNCTIONS
=#

# make objective function and objective function closure
begin
    function objective_function(x, p)
        geo = GeometricProperties_noadjoint(x, p)
        return dot(geo.A, geo.L)
    end

    OBJ = x -> objective_function(x, params)
end

# test objective function and gradient
@time o0, ∇o0 = withgradient(OBJ, x_init);

# make constraint function and closure
begin
    function constraint_function(x, p, dmax, fmax)
        res = solve_truss_noadjoint(x, p)

        vertical_displacements = res.U[3:3:end]

        flocal = [R * K * res.U[i] for (R, K, i) in zip(res.R, res.K, p.dofids)]
        axial_stresses = [abs(f[2]) / a for (f, a) in zip(flocal, res.A)]

        return [
            (-vertical_displacements .- dmax);
            (axial_stresses[i_stressed_elements] .- fmax);
        ]

    end

    CSTR = x -> constraint_function(x, params, dmax, fy)
end

# test constraint function and jacobian
@time c0, ∇c0 = withjacobian(CSTR, x_init);


# assert we start at a feasible point in the design space
@assert all(c0 .< 0) "CONSTRAINTS VIOLATED, INCREASE A0 TO START IN FEASIBLE POSITION"

#=
OPTIMIZATION
=#

alg = NLoptAlg(:LD_MMA)

opts = NLoptOptions(
    maxeval = 1000,
    maxtime = 300,
)

res = constrained_optimization(params, OBJ, CSTR, alg, opts)


# post process
begin
    @show res.obj_opt
    @show res.time
    @show n_iter = length(res)
    @show res.opt_stop_type

    x_opt = res.x_opt
    
    c_opt = CSTR(x_opt)

    c_baseline = [fill(dmax, model.nNodes); fill(fy, length(i_stressed_elements))]
    @show check_cstr(c_opt, c_baseline)
    feasibility = [check_cstr(c, c_baseline) for c in res.cstr_history]
    true_start, true_end, false_start, false_end = find_bands(feasibility)

    model_opt = updatemodel(params, x_opt)
    geo_opt = Geo(model_opt)
    dfac = Observable(0.)

    _nodes2 = @lift(Point3.(geo_opt.nodes .+ $dfac .* geo_opt.disp))
    _elements2 = @lift($_nodes2[geo_opt.indices_flat])
    
    lfac = Observable(4)
    lw = @lift(geo_opt.areas ./ geo_opt.max_area .* $lfac);

    tinc = range(0, res.time, n_iter)
    ts1 = tinc[true_start]
    ts2 = tinc[true_end]
    fs1 = tinc[false_start .- 1]
    fs2 = tinc[false_end]
end;

diameters = sqrt.(geo_opt.areas .* 4 ./ pi)
lw2 = diameters ./ maximum(diameters) .* 3

# visualize
begin
    fig = Figure(
        size = halfwidth()
    )
    ax = Axis3(
        fig[1,1],
        aspect = :data
    )

    hidedecorations!(ax); hidespines!(ax)

    linesegments!(
        _elements2,
        # color = geo_opt.forces,
        # colorrange = (-1, 1) .* geo_opt.max_abs_force .* .1,
        # colormap = pink2blue,
        linewidth = lw2,
    )

    # scatter!(
    #     _nodes2,
    #     color = :white,
    #     strokecolor = kjl_blue,
    #     strokewidth = 2
    # )

    axloss = Axis(
        fig[2,1],
        xlabel = "TIME [s]",
        ylabel = "MASS [tonnes]",
        aspect = nothing
    )

    style1!(axloss)

    vspan!(ts1, ts2, color = (kjl_blue, .15), label = "FEASIBLE")
    # vspan!(fs1, fs2, color = (kjl_pink, .15), label = "INFEASIBLE")

    lines!(
        tinc,
        res.obj_history .* 8, 
        color = :black,
        label = "OBJ"
        )

    leg = axislegend(axloss)

    fig
end

# jldsave(datadir * "07_28_2024_na_mma.jld2"; results = res)
# GHsave(model_opt, datadir * "07_28_2024_na_mma.json")