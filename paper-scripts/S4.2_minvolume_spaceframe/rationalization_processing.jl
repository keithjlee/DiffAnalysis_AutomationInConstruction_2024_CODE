using DiffAnalysis_AIC2024
set_theme!(aic)
include("init_problem.jl")

using CairoMakie; CairoMakie.activate!()
figdir = joinpath(@__DIR__, "figures/")

# read data
begin
    all_areas = jldopen(joinpath(@__DIR__, "data/07_28_2024_gb_mma.jld2"))["results"]
    three_areas = jldopen(joinpath(@__DIR__, "data/08_07_2024_gb_mma_3areas.jld2"))["results"]
end

model1 = all_areas.model_opt
model2 = three_areas.model_opt

geo1 = Geo(model1)
geo2 = Geo(model2)

tinc1 = range(0, all_areas.time, length(all_areas))
tinc2 = range(0, three_areas.time, length(three_areas))

begin
    fig = Figure(size = halfwidth())
    ax = Axis(
        fig[1,1],
        title = "OPTIMIZATION HISTORY",
        xlabel = "TIME [s]",
        ylabel = "VOLUME [m³]",
        aspect = 2.5
    )

    style1!(ax)
    lines!(tinc1, all_areas.obj_history, label = "GB-MMA")
    lines!(tinc2, three_areas.obj_history, label = "GB-MMA (3 AREAS)")

    leg = axislegend(ax, position = :rt)

    ax2 = Axis(
        fig[2,1],
        xlabel = "AREA [m²]",
        ylabel = "LENGTH [m]",
        aspect = 2.5,
        title = "ELEMENT PROPERTIES"
    )

    style1!(ax2)

    # vlines!(geo1.areas, color = (kjl_blue, .25))
    # vlines!(geo2.areas, color = (kjl_green, .25))

    ms = 7.5
    demand1 = scatter!(
        geo1.areas, geo1.lengths,
        markersize = ms,
        label = "GB-MMA",
        strokewidth = 0.5,
        strokecolor = :white
    )

    demand2 = scatter!(
        geo2.areas, geo2.lengths,
        markersize = ms,
        label = "GB-MMA (3 AREAS)",
        strokewidth = 0.5,
        strokecolor = :white
    )

    
    

    fig
end

save(figdir * "all_v_3areas.pdf", fig)