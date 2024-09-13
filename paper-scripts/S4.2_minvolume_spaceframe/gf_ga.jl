include("init_problem.jl")

# define objectives and constraints
begin
    function obj(x, p)
        geo = GeometricProperties(x, p)
        return dot(geo.L, geo.A)
    end

    OBJ = x -> obj(x, params)
    @time o0, do0 = withgradient(OBJ, x_init)

    function cstr(x, p, dmax, fmax)
        res = solve_truss(x, p)

        vertical_displacements = res.U[3:3:end]
        axial_stresses = axial_stress(res, p)

        return [
            (-vertical_displacements .- dmax);
            (axial_stresses[i_stressed_elements] .- fmax);
        ]
    end

    CSTR = x -> cstr(x, params, dmax, fy)
    @time c0, dc0 = withjacobian(CSTR, x_init)
end

@assert all(c0 .< 0)


#=
OPTIMIZATION
=#

# solver parameters
begin
    F = TraceFunction(OBJ)

    omodel = Nonconvex.Model(F)
    addvar!(omodel, params.lb, params.ub)
    add_ineq_constraint!(omodel, CSTR)

    Nonconvex.@load Metaheuristics

    alg = MetaheuristicsAlg(GA)
    opts = MetaheuristicsOptions(N = 100, f_calls_limit = 25_000)
end

begin
    t0 = time()
    res = Nonconvex.optimize(
        omodel,
        alg,
        x_init,
        options=opts
        );

    restime = time() - t0
end 

# post process
begin
    @show res.minimum
    @show restime
    @show n_iter = length(F.trace)

    x_opt = res.minimizer
    c_opt = CSTR(x_opt)
    
    c_baseline = [fill(dmax, model.nNodes); fill(fy, length(i_stressed_elements))]
    @show check_cstr(c_opt, c_baseline)
    

    x_history = getproperty.(F.trace, :input)
    cstr_history = [CSTR(x) for x in x_history]
    obj_history = getproperty.(F.trace, :output)

    feasibility = [check_cstr(c, c_baseline) for c in cstr_history]
    true_start, true_end, false_start, false_end = find_bands(feasibility)

    model_opt = updatemodel(params, x_opt)
    geo_opt = Geo(model_opt)

    _nodes2 = @lift(Point3.(geo_opt.nodes .+ $dfac .* geo_opt.disp))
    _elements2 = @lift($_nodes2[geo_opt.indices_flat])
    
    lfac = Observable(4)
    lw = @lift(geo_opt.areas ./ geo_opt.max_area .* $lfac);

    tinc = range(0, restime, n_iter)
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
    ax = LScene(
        fig[1,1],
        # aspect = :data
    )

    # hidespines!(ax); tickstoggle!(ax)
    # ylims!(-.5dy, 2.5dy)

    # linesegments!(
    #     _elements,
    #     # linestyle = :dash,
    #     linewidth = .5,
    #     color = (:black, .5)
    # )

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
    vspan!(fs1, fs2, color = (kjl_pink, .15), label = "INFEASIBLE")

    lines!(
        tinc,
        obj_history .* 8, 
        color = :black,
        label = "OBJ"
        )

    leg = axislegend(axloss)

    fig
end

# jldsave(
#     datadir * "07_28_2024_gf_ga.jld2";
#     time = restime,
#     params = params,
#     model_opt = model_opt,
#     x_opt = x_opt,
#     obj_opt = res.minimum,
#     cstr_opt = c_opt,
#     x_history = x_history,
#     obj_history = obj_history,
#     cstr_history = cstr_history,
#     alg = alg,
#     opts = opts,
#     opt_stop_type = :MAXEVAL_REACHED
# )

# GHsave(model_opt, datadir * "07_28_2024_gf_ga.json")