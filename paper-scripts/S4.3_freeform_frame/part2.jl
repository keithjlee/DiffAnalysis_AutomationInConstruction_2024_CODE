include("initialize2.jl")

#=
OPTIMIZATION
=#
# variables
begin
    strutvecs1 = Asap.local_x.(model.elements[:strut1]; unit = false)
    vector_parents = [VectorSpatialVariable(node, vec, 0., -0.5, 2.5) for (node, vec) in zip(model.nodes[:active1], strutvecs1)]

    strutvecs2 = Asap.local_x.(model.elements[:strut2]; unit = false)
    vector_children = [CoupledVariable(node, parent, vec) for (node, parent, vec) in zip(model.nodes[:active2], vector_parents, strutvecs2)]

    vars = FrameVariable[
        vector_parents;
        vector_children;
        [SpatialVariable(node, 0., -4., 0.1, :X) for node in model.nodes[:leftanchor]];
        [SpatialVariable(node, 0., -0.1, 4., :X) for node in model.nodes[:rightanchor]]
    ]
end

params = FrameOptParams(model, vars)
x0 = copy(params.values)

# objective
begin
    function obj(x, p)
        res = solve_frame(x, p)
        dot(res.U, p.P)
    end

    OBJ = x -> obj(x, params)
    @time o0, do0 = withgradient(OBJ, x0)
end


@time o0, do0 = withgradient(OBJ, x0)

# optimization params
begin
    alg = NLoptAlg(:LD_LBFGS)
    opts = NLoptOptions(
        maxeval = 1000,
        maxtime = 300,
        ftol_rel = 1e-5
    )
end

#solve
res = unconstrained_optimization(params, OBJ, alg, opts)

# post process
begin
    @show res.obj_opt
    @show res.time
    @show n_iter = length(res)
    @show res.opt_stop_type

    x_opt = res.x_opt
    c_opt = res.cstr_opt

    tinc = range(0, res.time, n_iter)

    gset = [Geo(updatemodel(params, x)) for x in res.x_history]

    i = Observable(1)
    g = @lift(gset[$i])
    p = @lift(Point3.($g.nodes))
    e = @lift($p[$g.indices_flat])
    lw = @lift($g.areas ./ maximum($g.areas) .* 3)
    l = @lift($g.lengths)

    m = @lift(abs.($g.Mz))
    loss = @lift(Point2.(tinc[1:$i], res.obj_history[1:$i]))
end

# visualize
begin
    fig = Figure(
        # size = halfwidth()
    )
    ax = LScene(fig[1,1], show_axis = false)

    linesegments!(
        e,
        linewidth = lw,
    )

    axloss = Axis(
        fig[2,1],
        xlabel = "TIME [s]",
        ylabel = "COMPLIANCE [kNm]",
        aspect = nothing
    )

    

    style1!(axloss)
    
    lines!(
        loss,
        color = :black,
        label = "OBJ"
        )

    on(i) do _
        autolimits!(axloss)
    end

    rowsize!(fig.layout, 2, Aspect(1, .3))

    fig
end

i[] = 1
for k = 1:length(res)
    i[] = k
    sleep(1e-2)
end

record(fig, "test.gif", 1:length(res); framerate = 20) do frame
    i[] = frame
end

# jldsave(datadir * "08_12_2024_spine2.jld2"; results = res)
# GHsave(res.model_opt, datadir * "08_12_2024_spine2.json")