include("init_problem.jl")

figdir = joinpath(@__DIR__, "figures/")
isdir(figdir) || mkdir(figdir)

using CairoMakie; CairoMakie.activate!()

begin
    fig = Figure(
        size = halfwidth(.5),
        backgroundcolor = :transparent
    )

    ax = Axis(
        fig[1,1],
        aspect = DataAspect()
    )

    hidedecorations!(ax)
    hidespines!(ax)

    ylims!(ax, -.5dy, 1.5dy)

    linesegments!(
        _elements,
        color = :black
    )

    scatter!(
        _nodes,
        strokewidth = 2,
        color = :white
    )

    fig
end

save(figdir * "init.pdf", fig)