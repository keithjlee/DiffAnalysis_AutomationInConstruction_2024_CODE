using DiffAnalysis_AIC2024
set_theme!(aic)
include("init_problem.jl")

using CairoMakie; CairoMakie.activate!()
figdir = joinpath(@__DIR__, "figures/")

# read data
begin
    gb_mma = jldopen(joinpath(@__DIR__, "data/07_25_2024_gb_mma.jld2"))["results"]
    fd_mma = jldopen(joinpath(@__DIR__, "data/07_25_2024_fd_mma.jld2"))["results"]
    na_mma = jldopen(joinpath(@__DIR__, "data/07_25_2024_na_mma.jld2"))["results"]

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
getproperty.(datasets, :time)

begin
    tinc_gb_mma = range(0, gb_mma.time, length(gb_mma))
    tinc_fd_mma = range(0, fd_mma.time, length(fd_mma))
    tinc_na_mma = range(0, na_mma.time, length(na_mma))
    tinc_gf_cobyla = range(0, gf_cobyla.time, length(gf_cobyla))
    tinc_gf_ga = range(0, gf_ga.time, length(gf_ga.x_history))
end

begin
    c_gb_mma = [all(c .< 0) for c in gb_mma.cstr_history]
    c_fd_mma = [all(c .< 0) for c in fd_mma.cstr_history]
    c_na_mma = [all(c .< 0) for c in na_mma.cstr_history]
    c_gf_cobyla = [all(c .< 0) for c in gf_cobyla.cstr_history]
    c_gf_ga = [all(c.< 0) for c in gf_ga.cstr_history]
end
#=
OPTIMIZATION HISTORY COMPARISON
=#

#overall optimization trace
begin
    lw = 2
    fig = Figure(size = halfwidth())

    ax = Axis(
        fig[1,1],
        ylabel = "VOLUME [m³]",
        aspect = nothing
    )

    style1!(ax)

    lines!(
        tinc_gf_cobyla[1:50:end], gf_cobyla.obj_history[1:50:end],
        linewidth = lw,
        label = "GF-COBYLA",
        linestyle = :dash,
        # color = :gray,
        color = c_gf_cobyla[1:50:end],
        colormap = cgrad([kjl_pink, :gray])
    )

    lines!(
        tinc_gf_ga[1:end], gf_ga.obj_history[1:end],
        linewidth = lw,
        label = "GF-GA",
        # color = :gray
        color = c_gf_ga[1:end],
        colormap = cgrad([kjl_pink, :gray])
    )

    lines!(
        tinc_fd_mma, fd_mma.obj_history,
        linewidth = lw,
        label = "FD-MMA",
        linestyle = :dash,
        # color = (kjl_blue, .75),
        color = c_fd_mma,
        colormap = cgrad([kjl_pink, kjl_blue]),
        colorrange = (0, 1)
    )

    lines!(
        tinc_na_mma, na_mma.obj_history,
        linewidth = lw,
        label = "NA-MMA",
        linestyle = :dot,
        # color = (kjl_blue, .75)
        color = c_na_mma,
        colormap = cgrad([kjl_pink, kjl_blue]),
        colorrange = (0, 1)
    )

    lines!(
        tinc_gb_mma, gb_mma.obj_history,
        linewidth = lw,
        label = "GB-MMA (OURS)",
        # color = kjl_blue
        color = c_gb_mma,
        colormap = cgrad([kjl_pink, kjl_blue]),
        colorrange = (0, 1)
    )

    xlim = 7.
    ylim = 4.75
    lines!(
        Rect2f(Point2(-1., 0.), Vec2(xlim, ylim)),
        color = :black,
        linewidth = 1
    )

    axislegend(ax)

    axzoom = Axis(
        fig[2,1],
        xlabel = "TIME [s]",
        ylabel = "VOLUME [m³]",
        aspect = nothing
    )

    xlims!(axzoom, nothing, xlim)
    ylims!(axzoom, nothing, ylim)

    style1!(axzoom)
    axzoom.topspinevisible = axzoom.rightspinevisible = true

    lines!(
        range(0, fd_mma.time, length(fd_mma))[1:5], fd_mma.obj_history[1:5],
        linewidth = lw,
        # label = "FD-MMA",
        linestyle = :dash,
        # color = (kjl_blue, .75)
        color = c_fd_mma[1:5],
        colormap = cgrad([kjl_pink, kjl_blue]),
        colorrange = (0, 1)
    )

    lines!(
        range(0, na_mma.time, length(na_mma))[1:60], na_mma.obj_history[1:60],
        linewidth = lw,
        # label = "NA-MMA",
        linestyle = :dot,
        # color = (kjl_blue, .75)
        color = c_na_mma[1:60],
        colormap = cgrad([kjl_pink, kjl_blue]),
        colorrange = (0, 1)
    )

    lines!(
        tinc_gf_cobyla[1:10:end], gf_cobyla.obj_history[1:10:end],
        linewidth = lw,
        label = "GF-COBYLA",
        linestyle = (:dash, :dense),
        color = c_gf_cobyla[1:10:end],
        colormap = cgrad([kjl_pink, :gray]),
        colorrange = (0, 1)
    )

    band!(
        tinc_gf_ga, gf_ga.obj_lower, gf_ga.obj_upper,
        color = (:gray, .25)
    )

    lines!(
        tinc_gf_ga, gf_ga.obj_history,
        linewidth = lw,
        label = "GF-GA",
        color = :gray
    )

    lines!(
        tinc_gb_mma, gb_mma.obj_history,
        linewidth = 2lw,
        label = "GB-MMA (OUR METHOD)",
        # color = kjl_blue
        color = c_gb_mma,
        colormap = cgrad([kjl_pink, kjl_blue]),
        colorrange = (0, 1)
    )

    xlim2 = 0.5
    ylim2 = 2.5
    lines!(
        Rect2f(Point2(0., 0.), Vec2(xlim2, ylim2)),
        color = :black,
        linewidth = 1
    )


    # axislegend(axzoom)

    # rowsize!(fig.layout, 1, Aspect(1, .5))

    fig
end

save(figdir * "optimization_trace_overall_v6.pdf", fig)

# zoom in either further for inset
begin
    fig = Figure(size = halfwidth())
    ax = Axis(
        fig[1,1],
        # aspect = 2
    )

    style1!(ax)
    ax.topspinevisible = ax.rightspinevisible = true

    xmin, xmax = 0., xlim2
    ymin, ymax = 0., ylim2

    hidedecorations!(ax)

    t_cobyla = range(0, gf_cobyla.time, length(gf_cobyla))
    t_ga = range(0, gf_ga.time, length(gf_ga.obj_history))
    t_mma = range(0, gb_mma.time, length(gb_mma))

    i_cobyla_stop = findfirst(t_cobyla .>= xmax)
    i_ga_stop = findfirst(t_ga .>= xmax)
    i_mma_stop = findfirst(t_mma .>= xmax)

    lines!(
        t_cobyla[1:2:i_cobyla_stop], gf_cobyla.obj_history[1:2:i_cobyla_stop],
        linewidth = lw,
        label = "GF-COBYLA",
        # color = :gray,
        linestyle = (:dash, :dense),
        color = c_gf_cobyla[1:2:i_cobyla_stop],
        colormap = cgrad([kjl_pink, :gray]),
        colorrange = (0, 1)
    )

    band!(
        tinc_gf_ga, gf_ga.obj_lower, gf_ga.obj_upper,
        color = (:gray, .25)
    )

    lines!(
        tinc_gf_ga, gf_ga.obj_history,
        linewidth = lw,
        label = "GF-GA",
        color = :gray
    )

    lines!(
        t_mma[1:i_mma_stop], gb_mma.obj_history[1:i_mma_stop],
        linewidth = 2lw,
        label = "GB-MMA (OUR METHOD)",
        color = kjl_blue
    )

    limits!(ax, xmin, xmax, ymin, ymax)

    fig
end

save(figdir * "optimization_trace_zoomed_v6.pdf", fig)

#=
FINAL SOLUTION GEOMETRIES
=#

begin
    geo_gb_mma = Geo(gb_mma.model_opt)
    geo_fd_mma = Geo(fd_mma.model_opt)
    geo_na_mma = Geo(na_mma.model_opt)

    geo_gf_cobyla = Geo(gf_cobyla.model_opt)
    geo_gf_ga = Geo(gf_ga.model_opt)
end

geos = [geo_gb_mma, geo_fd_mma, geo_na_mma, geo_gf_cobyla, geo_gf_ga]
names = ["gb_mma", "fd_mma", "na_mma", "gf_cobyla", "gf_ga"]

lfac = 4.
cfac = 1.

for (geo, name) in zip(geos, names)

    pts = Point2.(geo.nodes_xy)
    els = pts[geo.indices_flat]
    lw = geo.areas ./ geo.max_area .* lfac
    cr = geo.max_abs_force .* (-1, 1) .* cfac


    fig = Figure(backgroundcolor = :transparent)
    ax = Axis(
        fig[1,1], 
        aspect = DataAspect()
    )

    ylims!(-0.5, 2.5)
    hidedecorations!(ax); hidespines!(ax)

    # linesegments!(_elements, color = (:black, .5), linewidth = .5, linestyle = :dash)

    ls_elements = linesegments!(
        els,
        color = :black,
        linewidth = lw
    )

    sc_nodes = scatter!(
        pts,
        color = :white,
        strokewidth = 2,
        markersize = 10
    )

    save(figdir * name * "_sol1.pdf", fig)

    fig = Figure(
        backgroundcolor = :transparent
    )
    ax = Axis(
        fig[1,1], 
        aspect = DataAspect()
    )

    ylims!(-0.5, 2.5)
    hidedecorations!(ax); hidespines!(ax)

    ls_elements = linesegments!(
        els,
        color = geo.forces,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = lw,
    )

    sc_nodes = scatter!(
        pts,
        color = :white,
        strokewidth = 2,
        markersize = 10
    )

    save(figdir * name * "_sol2.pdf", fig)

    fig = Figure(
        backgroundcolor = :transparent
    )

    cb = Colorbar(
        fig[1,1],
        label = "AXIAL FORCE [kN]",
        vertical = false,
        flipaxis = false,
        colorrange = cr,
        colormap = pink2blue,
        tellheight = false
    )

    save(figdir * name * "_cb.pdf", fig)

end


#=
DEMAND SPACE
=#
GLMakie.activate!()
names = ["GB-MMA", "FD-MMA", "NA-MMA", "GF-COBYLA", "GF-GA"]
begin
    fig = Figure(size = halfwidth())

    ax = Axis(fig[1,1],
        xlabel = "AREA [m²]",
        ylabel = "LENGTH [m]"
    )

    style1!(ax)

    for (geo, name) in zip(geos, names)
        scatter!(geo.areas, geo.lengths, label = name, markersize = 15)
    end

    axislegend(ax)
    fig

end

#=
BI-OBJECTIVE
=#

times = getproperty.(datasets, :time)
volumes = getproperty.(datasets, :obj_opt)

tv_points = Point2.(times, volumes)
validity = [all(data.cstr_opt .< 0) for data in datasets]
colors = [kjl_blue, RGBf(0.,0.,0.), kjl_pink, kjl_pink, RGBf(0., 0., 0.)]

begin
    fig = Figure(size = halfwidth())

    ax = Axis(
        fig[1,1],
        xlabel = "TIME [s]",
        ylabel = "FINAL VOLUME [m³]",
        # aspect = nothing
        aspect = 2
    )

    style1!(ax)

    scatter!(tv_points, strokecolor = colors, color = :white, markersize = 15, strokewidth = 2)


    fig
end

# save(figdir * "time_v_sol.pdf", fig)