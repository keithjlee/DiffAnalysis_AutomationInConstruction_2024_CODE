include("init.jl")


geo = Geo(model)
pts = Point2.(geo.nodes_xy)

els_bottom = pts[vcat(geo.indices[i_bottom_elements]...)]
els_top = pts[vcat(geo.indices[i_top_elements]...)]
els_bottomweb = pts[vcat(geo.indices[i_bottomweb_elements]...)]
els_topweb = pts[vcat(geo.indices[i_topweb_elements]...)]
els_strut = pts[vcat(geo.indices[i_strut_elements]...)]
els_support = pts[vcat(geo.indices[i_support_elements]...)]

els_set = [els_bottom, els_top, els_bottomweb, els_topweb, els_strut, els_support]


begin
    fig = Figure(backgroundcolor = :transparent, size = halfwidth())
    ax = Axis(fig[1,1], aspect = DataAspect())
    ax.titlesize = 28

    hidedecorations!(ax); hidespines!(ax)

    # ax.titlevisible = true
    ylims!(ax, -dy, 2dy + 10 + support_vertical_offset)

    linesegments!(static_deck_points, color = :gray, linewidth = 1)

    ls_elements = [linesegments!(
        e,
        color = :black,
        linewidth = 1
    ) for e in els_set]

    scatter!(pts, color = :white, markersize = 3, strokewidth = .5)

    fig
end

save(figdir * "init.pdf", fig)

begin
    fig = Figure(backgroundcolor = :transparent, size = halfwidth())
    ax = Axis(fig[1,1], aspect = DataAspect())
    ax.titlesize = 28

    hidedecorations!(ax); hidespines!(ax)

    # ax.titlevisible = true
    ylims!(ax, -dy, 2dy + 10 + support_vertical_offset)

    linesegments!(static_deck_points, color = :gray, linewidth = 1)

    ls_elements = [linesegments!(
        e,
        linewidth = 1
    ) for e in els_set]

    scatter!(pts, color = :white, markersize = 3, strokewidth = .5)

    fig
end

save(figdir * "init_colored.pdf", fig)