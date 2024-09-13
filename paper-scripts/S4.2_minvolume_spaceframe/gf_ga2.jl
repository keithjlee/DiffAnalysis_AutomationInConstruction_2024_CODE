include("init_problem.jl")
using Metaheuristics
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

function OBJ_CSTR(x)
    return OBJ(x), CSTR(x), [0.]
end

OBJ_CSTR(x_init)

#=
OPTIMIZATION
=#
#Metaheuristics
begin
    bounds = BoxConstrainedSpace(lb = params.lb, ub = params.ub)
    opts = Metaheuristics.Options(time_limit = 300.0, store_convergence = true, f_tol = 1e-6, f_calls_limit = 200_000, iterations = 2000)
    alg = Metaheuristics.GA(N = 100, options = opts, mutation = PolynomialMutation(;bounds))
end

begin
    t0 = time()
    res = Metaheuristics.optimize(OBJ_CSTR, bounds, alg)
    restime = time() - t0
end 


# post process
begin
    @show res.best_sol.f
    @show restime = res.overall_time
    @show n_iter = res.iteration
    x_opt = res.best_sol.x
    c_opt= CSTR(x_opt)
    c_baseline = [fill(dmax, model.nNodes); fill(fy, model.nElements)]

    tinc = range(0, restime, n_iter)
end;

#histories
begin
    x_history = Vector{Vector{Float64}}()
    obj_history = Vector{Float64}()
    obj_upper = Vector{Float64}()
    obj_lower = Vector{Float64}()
    cstr_history = Vector{Vector{Float64}}()
    feasibility = Vector{Bool}()

    for iter in res.convergence

        push!(x_history, iter.best_sol.x)
        push!(obj_history, iter.best_sol.f)
        push!(cstr_history, iter.best_sol.g)
        push!(feasibility, iter.best_sol.is_feasible)

        #bounds
        fvals = getproperty.(iter.population, :f)
        push!(obj_lower, minimum(fvals))
        push!(obj_upper, maximum(fvals))

    end
end

begin
    model_opt = updatemodel(params, x_opt)
    geo_opt = Geo(model_opt)

    _nodes2 = @lift(Point3.(geo_opt.nodes .+ $dfac .* geo_opt.disp))
    _elements2 = @lift($_nodes2[geo_opt.indices_flat])
    
    lfac = Observable(4)
    lw = @lift(geo_opt.areas ./ geo_opt.max_area .* $lfac);
end

# visualize

begin
    fig = Figure(
        size = halfwidth()
    )
    ax = Axis3(
        fig[1,1],
        aspect = :data
    )

    hidespines!(ax); tickstoggle!(ax)

    # linesegments!(
    #     _elements,
    #     linestyle = :dash,
    #     linewidth = .5,
    #     color = :black
    # )

    linesegments!(
        _elements2,
        # color = geo_opt.forces,
        # colorrange = (-1, 1) .* geo_opt.max_abs_force .* .25,
        # colormap = pink2blue,
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

    band!(tinc, obj_lower .* 8, obj_upper .* 8)

    lines!(
        tinc,
        obj_history .* 8, 
        color = :black,
        label = "OBJ"
        )


    fig
end

# jldsave(
#     datadir * "07_29_2024_gf_ga.jld2";
#     time = restime,
#     params = params,
#     model_opt = model_opt,
#     x_opt = x_opt,
#     obj_opt = res.best_sol.f,
#     cstr_opt = c_opt,
#     x_history = x_history,
#     obj_history = obj_history,
#     obj_lower = obj_lower,
#     obj_upper = obj_upper,
#     cstr_history = cstr_history,
#     n_f_calls = res.f_calls,
#     alg = alg,
#     opts = opts,
#     opt_stop_type = :MAXEVAL_REACHED
# )

# GHsave(model_opt, datadir * "07_29_2024_gf_ga.json")