using DiffAnalysis_AIC2024
# optimization trace for shape optimization
figdir = joinpath(@__DIR__, "figures/")
isdir(figdir) || mkdir(figdir)

shapeopt_data = jldopen(joinpath(@__DIR__, "data/08_12_2024_spine2.jld2"))["results"]

tinc = range(0, shapeopt_data.time, length(shapeopt_data))
loss = shapeopt_data.obj_history

using CairoMakie; CairoMakie.activate!()
begin
    fig = Figure(size = halfwidth(0.5))

    ax = Axis(
        fig[1,1],
        aspect = nothing,
        xlabel = "TIME [s]",
        ylabel = "COMPLIANCE [kNm]"
    )

    style1!(ax)

    lines!(tinc, loss, color = :black)

    fig
end

save(figdir * "shapeopt.pdf", fig)

# size opt
sizeopt_data = jldopen(joinpath(@__DIR__, "data/08_12_2024_final2.jld2"))

tinc_size = sizeopt_data["time_history"]
loss_size = sizeopt_data["obj_history"]

begin
    fig = Figure(size = halfwidth(0.5))

    ax = Axis(
        fig[1,1],
        aspect = nothing,
        xlabel = "TIME [s]",
        ylabel = "VOLUME [m³]"
    )

    style1!(ax)

    lines!(tinc_size, loss_size, color = :black)

    fig
end

save(figdir * "sizeopt.pdf", fig)

dframe = sizeopt_data["dframe"]
tframe = (dframe * (1 - sizeopt_data["aframe"])) / 2

dstrut = sizeopt_data["dstrut"]
tstrut = (dstrut * (1 - sizeopt_data["astrut"])) / 2

dspine = sizeopt_data["dspine"]
tspine = (dspine * (1 - sizeopt_data["aspine"])) / 2

dtie = sizeopt_data["dtie"]
ttie = (dtie * (1 - sizeopt_data["atie"])) / 2

