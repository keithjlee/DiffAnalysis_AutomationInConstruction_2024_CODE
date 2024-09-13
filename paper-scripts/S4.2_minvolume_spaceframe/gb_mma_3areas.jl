include("init_problem_3areas.jl")

begin
    fig = Figure()
    ax = Axis3(
        fig[1,1],
        aspect = :data,
    )

    linesegments!(_elements, linewidth = 1, color = :black)

    fig
end

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

@assert all(c0 .< 0)

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
end;

# visualize
begin
    fig = Figure(
        size = halfwidth()
    )
    ax = LScene(
        fig[1,1],
    )

    linesegments!(
        _elements2,
        linewidth = lw,
    )

    axloss = Axis(
        fig[2,1],
        xlabel = "TIME [s]",
        ylabel = "MASS [tonnes]",
        aspect = nothing
    )

    style1!(axloss)

    lines!(
        tinc,
        res.obj_history .* 8, 
        color = :black,
        label = "OBJ"
        )

    leg = axislegend(axloss)

    fig
end

# jldsave(datadir * "08_07_2024_gb_mma_3areas.jld2"; results = res)
# GHsave(model_opt, datadir * "08_07_2024_gb_mma_3areas.json")