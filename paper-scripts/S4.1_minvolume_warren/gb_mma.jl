include("init_problem.jl")

#=
OPTIMIZATION: FUNCTIONS
=#

# make objective function and objective function closure
begin
    function objective_function(x, p)
        geo = GeometricProperties(x, p)
        return dot(geo.A, geo.L)
    end

    OBJ = x -> objective_function(x, params)
end

# test objective function and gradient
@time o0, do0 = withgradient(OBJ, x_init)

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

    end

    CSTR = x -> constraint_function(x, params, dmax, fy)
end

# test constraint function and jacobian
@time c0, dc0 = withjacobian(CSTR, x_init)

# assert we start at a feasible point in the design space
@assert all(c0 .< 0) "CONSTRAINTS VIOLATED, INCREASE A0 TO START IN FEASIBLE POSITION"

#=
OPTIMIZATION
=#

# solver parameters
begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(
        maxeval = 2000,
        maxtime = 12
    )
end

# solve
res = constrained_optimization(params, OBJ, CSTR, alg, opts);

# post process
begin
    @show res.obj_opt
    @show res.time
    @show n_iter = length(res)
    @show res.opt_stop_type

    x_opt = res.x_opt
    c_opt = CSTR(x_opt)

    c_baseline = [fill(dmax, model.nNodes); fill(fy, model.nElements)]
    @show check_cstr(c_opt, c_baseline)
    feasibility = [check_cstr(c, c_baseline) for c in res.cstr_history]
    true_start, true_end, false_start, false_end = find_bands(feasibility)

    model_opt = updatemodel(params, x_opt)
    geo_opt = Geo(model_opt)

    _nodes2 = @lift(Point2.(geo_opt.nodes_xy .+ $dfac .* geo_opt.disp_xy))
    _elements2 = @lift($_nodes2[geo_opt.indices_flat])
    
    lfac = Observable(4)
    lw = @lift(geo_opt.areas ./ geo_opt.max_area .* $lfac);

    tinc = range(0, res.time, n_iter)
    ts1 = tinc[true_start]
    ts2 = tinc[true_end]
    fs1 = tinc[false_start .- 1]
    fs2 = tinc[false_end]
end;

# visualize
begin
    fig = Figure(
        size = halfwidth()
    )
    ax = Axis(
        fig[1,1],
        aspect = DataAspect()
    )

    hidespines!(ax); tickstoggle!(ax)
    ylims!(-.5dy, 2.5dy)

    linesegments!(
        _elements,
        linestyle = :dash,
        linewidth = .5,
        color = :black
    )

    linesegments!(
        _elements2,
        linewidth = lw,
    )

    scatter!(
        _nodes2,
        color = :white,
        strokecolor = kjl_blue,
        strokewidth = 2
    )

    axloss = Axis(
        fig[2,1],
        xlabel = "TIME [s]",
        ylabel = "MASS [tonnes]",
        aspect = nothing
    )

    style1!(axloss)

    vspan!(ts1, ts2, color = (kjl_blue, .15), label = "FEASIBLE")
    vspan!(fs1, fs2, color = (kjl_pink, .15), label = "INFEASIBLE")

    lines!(
        tinc,
        res.obj_history .* 8, 
        color = :black,
        label = "OBJ"
        )

    leg = axislegend(axloss)

    fig
end

# save data
# jldsave(datadir * "07_25_2024_gb_mma.jld2"; results = res)
# GHsave(model_opt, datadir * "07_25_2024_gb_mma.json")