"""
    discretize(n::Integer; colormap = :tempo)

Discretize a color gradient into `n` discrete colours. Returns a vector of colors that can be indexed from 1 to `n`.
"""
function discretize(n::Integer; colormap = :tempo)
    return [cgrad(colormap, [0., 1.])[z] for z ∈ range(0., 1., length = n)]
end

"""
    labelscale!(axis::Axis, factor::Union{Float64, Int64})

Scale the font size of the title and labels by `factor`
"""
function labelscale!(axis::Axis, factor::Union{Float64, Int64})
    axis.xlabelsize = labelFontSize * factor
    axis.ylabelsize = labelFontSize * factor
    axis.xticklabelsize = tickFontSize * factor
    axis.yticklabelsize = tickFontSize * factor
    axis.titlesize = titleFontSize * factor
end

"""
    labelscale!(axis::Axis3, factor::Union{Float64, Int64})

Scale the font size of the title and labels by `factor`
"""
function labelscale!(axis::Axis3, factor::Union{Float64, Int64})
    axis.titlesize = titleFontSize * factor
    axis.xlabelsize = labelFontSize * factor
    axis.ylabelsize = labelFontSize * factor
    axis.zlabelsize = labelFontSize * factor
    axis.xticklabelsize = tickFontSize * factor
    axis.yticklabelsize = tickFontSize * factor
    axis.zticklabelsize = tickFontSize * factor
end

"""
    resetlabelscale!(axis::Axis)

Reset the font size to default:
- tick labels = 18
- axis labels = 20
- title = 20
"""
function resetlabelscale!(axis::Axis)
    axis.xlabelsize = labelFontSize
    axis.ylabelsize = labelFontSize
    axis.xticklabelsize = tickFontSize
    axis.yticklabelsize = tickFontSize
    axis.titlesize = titleFontSize
end

"""
    resetlabelscale!(axis::Axis3)

Reset the font size to default:
- tick labels = 18
- axis labels = 20
- title = 20
"""
function resetlabelscale!(axis::Axis3)
    axis.titlesize = titleFontSize
    axis.xlabelsize = labelFontSize
    axis.ylabelsize = labelFontSize
    axis.zlabelsize = labelFontSize
    axis.xticklabelsize = tickFontSize
    axis.yticklabelsize = tickFontSize
    axis.zticklabelsize = tickFontSize
end


"""
    gridtoggle!(axis::Axis)

Toggle the visibility of the grid
"""
function gridtoggle!(axis::Axis)
    axis.xgridvisible = !axis.xgridvisible[]
    axis.ygridvisible = !axis.ygridvisible[]
end

"""
    gridtoggle!(axis::Axis3)

Toggle the visibility of the grid
"""
function gridtoggle!(axis::Axis3)
    axis.xgridvisible = !axis.xgridvisible[]
    axis.ygridvisible = !axis.ygridvisible[]
    axis.zgridvisible = !axis.zgridvisible[]
end

"""
    simplifyspines!(axis::Axis3)

Simplify spines of an Axis3 to have one x/y/z spine
"""
function simplifyspines!(axis::Axis3)

    if axis.xspinecolor_2 != :transparent
        axis.xspinecolor_2 = :transparent
        axis.xspinecolor_3 = :transparent

        axis.yspinecolor_2 = :transparent
        axis.yspinecolor_3 = :transparent

        axis.zspinecolor_2 = :transparent
        axis.zspinecolor_3 = :transparent
    else
        axis.xspinecolor_2 = axis.xspinecolor_1[]
        axis.xspinecolor_3 = axis.xspinecolor_1[]

        axis.yspinecolor_2 = axis.xspinecolor_1[]
        axis.yspinecolor_3 = axis.xspinecolor_1[]

        axis.zspinecolor_2 = axis.xspinecolor_1[]
        axis.zspinecolor_3 = axis.xspinecolor_1[]
    end
    
end

"""
    mirrorticks!(axis::Axis)

Add mirrored ticks on the top and right spines
"""
function mirrorticks!(axis::Axis)
    axis.xticksmirrored = !axis.xticksmirrored[]
    axis.yticksmirrored = !axis.yticksmirrored[]
end

"""
    alignticks!(axis::Axis, value::Integer)

Position of ticks: 0 for outside, 1 for inside
"""
function alignticks!(axis::Axis, value::Integer)
    @assert value == 0 || value == 1 "Value must be 0 or 1"

    axis.xtickalign = value
    axis.ytickalign = value
end

"""
    tickstoggle!(axis::Union{Axis, Axis3})

Toggle the visibility of ticks
"""
function tickstoggle!(axis::Union{Axis, Axis3})
    axis.xticksvisible = !axis.xticksvisible[]
    axis.yticksvisible = !axis.yticksvisible[]

    if typeof(axis) == Axis3
        axis.zticksvisible = !axis.zticksvisible[]
    end
end

"""
    Get the size of a figure in pts
"""
Base.size(fig::Figure) = fig.scene.viewport.val.widths

"""
    fixlimits!(ax::Axis)

Fix axis limits to the current state
"""
function fixlimits!(ax::Axis)

    lx, ly = copy(ax.finallimits[].origin)
    ux, uy = copy(ax.finallimits[].widths)

    ax.limits = (lx, ux, ly, uy)

end

"""
    fixlimits!(ax::Axis3)
    
Fix axis limits to the current state
"""
function fixlimits!(ax::Axis3)

    lx, ly, lz = copy(ax.finallimits[].origin)
    ux, uy, uz = copy(ax.finallimits[].widths)

    ax.limits = (lx, ux, ly, uy, lz, uz)

end


function style1!(axis::Axis; bgcolor = :lightgray, bgalpha = 0.15, grid = true, gridcolor = :white)

    axis.backgroundcolor = (bgcolor, bgalpha)
    axis.xticksvisible = false
    axis.yticksvisible = false
    axis.rightspinevisible = false
    axis.topspinevisible = false

    if grid
        axis.xgridvisible = axis.ygridvisible = true
        axis.xgridcolor = axis.ygridcolor = gridcolor
    end

end

function style1!(axis::Axis3; bgcolor = :lightgray, bgalpha = 0.15, grid = false, gridcolor = :white)

    axis.xzpanelcolor = axis.yzpanelcolor = axis.xypanelcolor = (bgcolor, bgalpha)
    axis.xticksvisible = false
    axis.yticksvisible = false
    axis.zticksvisible = false

    if grid
        axis.xgridvisible = axis.ygridvisible = axis.zgridvisible = true
        axis.xgridcolor = axis.ygridcolor = axis.zgridcolor = gridcolor
    end

end