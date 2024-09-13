using DiffAnalysis_AIC2024
set_theme!(aic)
include("init_problem.jl")

using CairoMakie; CairoMakie.activate!()
figdir = joinpath(@__DIR__, "figures/")

# read data
begin
    gb_mma = jldopen(joinpath(@__DIR__, "data/07_28_2024_gb_mma.jld2"))["results"]
    fd_mma = jldopen(joinpath(@__DIR__, "data/07_29_2024_fd_mma.jld2"))["results"]
    na_mma = jldopen(joinpath(@__DIR__, "data/07_28_2024_na_mma.jld2"))["results"]

    gf_cobyla = jldopen(joinpath(@__DIR__, "data/07_29_2024_gf_cobyla.jld2"))["results"]
    gf_ga_data = jldopen(joinpath(@__DIR__, "data/07_29_2024_gf_ga.jld2"))
    gf_ga = (
        time = gf_ga_data["time"],
        params = gf_ga_data["params"],
        model_opt = gf_ga_data["model_opt"],
        x_opt = gf_ga_data["x_opt"],
        obj_opt = gf_ga_data["obj_opt"],
        cstr_opt = gf_ga_data["cstr_opt"],
        x_history = gf_ga_data["x_history"],
        obj_history = gf_ga_data["obj_history"],
        cstr_history = gf_ga_data["cstr_history"],
        alg = gf_ga_data["alg"],
        opts = gf_ga_data["opts"],
        opt_stop_type = gf_ga_data["opt_stop_type"],
        obj_lower = gf_ga_data["obj_lower"],
        obj_upper = gf_ga_data["obj_upper"]
    )
end

datasets = [gb_mma, fd_mma, na_mma, gf_cobyla, gf_ga]
setnames = ["GB-MMA", "FD-MMA", "NA-MMA", "GF-COBYLA", "GF-GA"]
getproperty.(datasets, :time)

[setnames getproperty.(datasets, :time)]

#=
OPTIMIZATION HISTORY COMPARISON
=#

#time increments
begin
    tinc_gf_cobyla = range(0, gf_cobyla.time, length(gf_cobyla))
    tinc_gf_ga = range(0, gf_ga.time, length(gf_ga.obj_history))

    tinc_fd_mma = range(0, fd_mma.time, length(fd_mma))
    tinc_na_mma = range(0, na_mma.time, length(na_mma))
    tinc_gb_mma = range(0, gb_mma.time, length(gb_mma))
end

c_gb_mma = [all(c .< 0) for c in gb_mma.cstr_history]
c_fd_mma = [all(c .< 0) for c in fd_mma.cstr_history]
c_na_mma = [all(c .< 0) for c in na_mma.cstr_history]
c_gf_cobyla = [all(c .< 0) for c in gf_cobyla.cstr_history]
c_gf_ga = [all(c.< 0) for c in gf_ga.cstr_history]

#overall optimization trace
# GLMakie.activate!()
begin
    lw = 2
    fig = Figure(size = halfwidth())

    ax = Axis(
        fig[1,1],
        xlabel = "TIME [s]",
        ylabel = "VOLUME [m³]",
        aspect = 2
    )

    style1!(ax)

    lines!(
        tinc_gf_cobyla, gf_cobyla.obj_history,
        linewidth = lw,
        label = "GF-COBYLA",
        color = :gray,
        linestyle = :dash
    )

    band!(
        tinc_gf_ga, gf_ga.obj_lower, gf_ga.obj_upper,
        color = (:gray, 0.25)
    )

    lines!(
        tinc_gf_ga, gf_ga.obj_history,
        linewidth = lw,
        label = "GF-GA",
        color = :gray
    )

    lines!(
        tinc_fd_mma, fd_mma.obj_history,
        linewidth = lw,
        label = "FD-MMA",
        linestyle = :dash,
        color = (kjl_blue, .75)
    )

    lines!(
        tinc_na_mma, na_mma.obj_history,
        linewidth = lw,
        label = "NA-MMA",
        linestyle = :dot,
        color = (kjl_blue, .75)
    )

    lines!(
        tinc_gb_mma, gb_mma.obj_history,
        linewidth = lw,
        label = "GB-MMA (OURS)",
        color = kjl_blue
    )

    axislegend(ax)

    fig
end

save(figdir * "obj_traces_wide.pdf", fig)

#=
analyze crossover solution between GA and GB-MMA
=#
begin
    obj_distmatrix = [abs(o1 - o2) for o1 in gb_mma.obj_history, o2 in gf_ga.obj_history]
    i_crossover = argmin(obj_distmatrix)
    i_gb_mma_crossover = i_crossover[1]
    i_gf_ga_crossover = i_crossover[2]

    x_gb_mma_crossover = gb_mma.x_history[i_gb_mma_crossover]
    x_gf_ga_crossover = gf_ga.x_history[i_gf_ga_crossover]

    model_gb_mma_crossover = updatemodel(gb_mma.params, x_gb_mma_crossover)
    model_gf_ga_crossover = updatemodel(gf_ga.params, x_gf_ga_crossover)

    geo_gb_mma_crossover = Geo(model_gb_mma_crossover)
    p_gb_mma_crossover = Point3.(geo_gb_mma_crossover.nodes)
    e_gb_mma_crossover = p_gb_mma_crossover[geo_gb_mma_crossover.indices_flat]
    lw_gb_mma_crossover = geo_gb_mma_crossover.areas ./ maximum(geo_gb_mma_crossover.areas) .* 4

    geo_gf_ga_crossover = Geo(model_gf_ga_crossover)
    p_gf_ga_crossover = Point3.(geo_gf_ga_crossover.nodes)
    e_gf_ga_crossover = p_gf_ga_crossover[geo_gf_ga_crossover.indices_flat]
    lw_gf_ga_crossover = geo_gf_ga_crossover.areas ./ maximum(geo_gf_ga_crossover.areas) .* 4
end

begin
    fig = Figure()
    ax = Axis3(
        fig[1,1],
        aspect = :data
    )

    linesegments!(
        e_gb_mma_crossover,
        linewidth = lw_gb_mma_crossover
    )


    ax = Axis3(
        fig[1,2],
        aspect = :data
    )

    linesegments!(
        e_gf_ga_crossover,
        linewidth = lw_gf_ga_crossover
    )

    fig

end

# GHsave(model_gb_mma_crossover, datadir * "07_29_2024_gb_mma_crossover.json")
# GHsave(model_gf_ga_crossover, datadir * "07_29_2024_gf_ga_crossover.json")

#=
BI-OBJECTIVE
=#

times = getproperty.(datasets, :time)
volumes = getproperty.(datasets, :obj_opt)

tv_points = Point2.(times, volumes)
validity = [all(data.cstr_opt .< 0) for data in datasets]
colors = [kjl_blue, RGBf(0.,0.,0.), RGBf(0.,0.,0.), kjl_pink, RGBf(0., 0., 0.)]

begin
    fig = Figure(size = halfwidth())

    ax = Axis(
        fig[1,1],
        xlabel = "TIME [s]",
        ylabel = "FINAL VOLUME [m³]",
        # aspect = nothing
    )

    style1!(ax)

    scatter!(tv_points, strokecolor = colors, color = :white, markersize = 15, strokewidth = 2)


    fig
end


#=

=#