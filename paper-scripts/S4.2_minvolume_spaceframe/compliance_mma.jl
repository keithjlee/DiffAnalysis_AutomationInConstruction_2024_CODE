#=
Investigate the result when we consider a minimum compliance problem instead
=#

include("init_problem.jl")

#=
Make special parmeters: e.g. start with very small areas
=#
begin
    vars = [SpatialVariable(node, 0. -.5*dz, 2dz, :Z) for node in model.nodes[:top]]
    vars = Vector{TrussVariable}()
    # independent variable diagonals
    for i in diag(gen.itop)
        if model.nodes[i].id == :support
            continue
        end

        push!(vars, SpatialVariable(model.nodes[i], 0., zmintop, zmaxtop, :Z))
    end

    # top variables coupled
    for i = 2:size(gen.itop, 1)
        for j = 1:i-1

            iparent = gen.itop[i,j]
            ichild = gen.itop[j,i]

            if model.nodes[iparent].id == :support || model.nodes[ichild].id == :support
                continue
            end

            parent = SpatialVariable(model.nodes[iparent], 0., zmintop, zmaxtop, :Z)
            child = CoupledVariable(model.nodes[ichild], parent)

            push!(vars, parent, child)

        end
    end

    for element in model.elements
        push!(vars, AreaVariable(element, .25Amax, Amin, Amax))
    end

    params = TrussOptParams(model, vars)
    x_init = params.values
end

res0 = solve_truss(x_init, params)
v0 = dot(res0.A, res0.L)

c0 = dot(res0.U, params.P)
std0 = std(abs.(res0.U[3:3:end]))

#=
OPTIMIZATION: FUNCTIONS
=#

# make objective function and objective function closure
begin
    function objective_function(x, p)
        res = solve_truss(x, p)
        std(abs.(res.U[3:3:end])) / std0 + dot(res.U, p.P) / c0
    end

    function constraint_function(x, p, vmax)
        geo = GeometricProperties(x, p)
        dot(geo.A, geo.L) - vmax
    end

    OBJ = x -> objective_function(x, params)
    CSTR = x -> constraint_function(x, params, .25v0)
end

# test objective function and gradient
@time o0, do0 = withgradient(OBJ, x_init);
@time c0, dc0 = withgradient(CSTR, x_init)

#=
OPTIMIZATION
=#

# solver parameters
begin
    alg = NLoptAlg(:LD_MMA)
    opts = NLoptOptions(
        maxeval = 2000,
        maxtime = 120,
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

    model_opt = updatemodel(params, x_opt)
    geo_opt = Geo(model_opt)

    _nodes2 = @lift(Point3.(geo_opt.nodes .+ $dfac .* geo_opt.disp))
    _elements2 = @lift($_nodes2[geo_opt.indices_flat])
    
    lfac = Observable(4)
    lw = @lift(geo_opt.areas ./ geo_opt.max_area .* $lfac);

    tinc = range(0, res.time, n_iter)
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
        ylabel = "COMPLIANCE [kNm]",
        aspect = nothing
    )

    style1!(axloss)
    
    lines!(
        tinc,
        res.obj_history .* 8, 
        color = :black,
        label = "OBJ"
        )


    fig
end

v2 = Asap.volume(model_opt)
stresses = geo_opt.forces ./ geo_opt.areas
extrema(stresses)

zdisps = getindex.(geo_opt.disp, 3)
extrema(zdisps)

@show extrema(stresses)
@show extrema(zdisps)