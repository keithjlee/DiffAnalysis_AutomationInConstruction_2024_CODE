include("TrussIndexer.jl")
include("NetworkIndexer.jl")
include("FrameIndexer.jl")

function populate!(indexer::AbstractIndexer, var::SpatialVariable{SpatialX})

    push!(indexer.iX, var.i)
    push!(indexer.iXg, var.iglobal)
    push!(indexer.fX, 1.0)

    indexer.activeX = true
end

function populate!(indexer::AbstractIndexer, var::SpatialVariable{SpatialY})

    push!(indexer.iY, var.i)
    push!(indexer.iYg, var.iglobal)
    push!(indexer.fY, 1.0)

    indexer.activeY = true
end

function populate!(indexer::AbstractIndexer, var::SpatialVariable{SpatialZ})

    push!(indexer.iZ, var.i)
    push!(indexer.iZg, var.iglobal)
    push!(indexer.fZ, 1.0)

    indexer.activeZ = true
end

function populate!(indexer::AbstractIndexer, var::VectorSpatialVariable)

    x, y, z = var.vec

    if x != 0.0
        push!(indexer.iX, var.i)
        push!(indexer.iXg, var.iglobal)
        push!(indexer.fX, x)

        indexer.activeX = true
    end

    if y != 0.0
        push!(indexer.iY, var.i)
        push!(indexer.iYg, var.iglobal)
        push!(indexer.fY, y)

        indexer.activeY = true
    end

    if z != 0.0
        push!(indexer.iZ, var.i)
        push!(indexer.iZg, var.iglobal)
        push!(indexer.fZ, z)

        indexer.activeZ = true
    end

end

function populate!(indexer::AbstractIndexer, var::CoupledVariable{VectorSpatialVariable})

    x, y, z = var.factor

    if x != 0.0
        push!(indexer.iX, var.i)
        push!(indexer.iXg, var.iglobal)
        push!(indexer.fX, x)

        indexer.activeX = true
    end

    if y != 0.0
        push!(indexer.iY, var.i)
        push!(indexer.iYg, var.iglobal)
        push!(indexer.fY, y)

        indexer.activeY = true
    end

    if z != 0.0
        push!(indexer.iZ, var.i)
        push!(indexer.iZg, var.iglobal)
        push!(indexer.fZ, z)

        indexer.activeZ = true
    end
end

function populate!(indexer::Union{TrussOptIndexer, FrameOptIndexer}, var::AreaVariable)

    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, 1.0)

    indexer.activeA = true

end

function populate!(indexer::AbstractIndexer, var::CoupledVariable{SpatialVariable{SpatialX}})

    push!(indexer.iX, var.i)
    push!(indexer.iXg, var.iglobal)
    push!(indexer.fX, var.factor)

end

function populate!(indexer::AbstractIndexer, var::CoupledVariable{SpatialVariable{SpatialY}})

    push!(indexer.iY, var.i)
    push!(indexer.iYg, var.iglobal)
    push!(indexer.fY, var.factor)

end

function populate!(indexer::AbstractIndexer, var::CoupledVariable{SpatialVariable{SpatialZ}})

    push!(indexer.iZ, var.i)
    push!(indexer.iZg, var.iglobal)
    push!(indexer.fZ, var.factor)
    
end

function populate!(indexer::Union{TrussOptIndexer, FrameOptIndexer}, var::CoupledVariable{AreaVariable})

    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionA}})

    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionIx}})

    push!(indexer.iIx, var.i)
    push!(indexer.iIxg, var.iglobal)
    push!(indexer.fIx, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionIy}})

    push!(indexer.iIy, var.i)
    push!(indexer.iIyg, var.iglobal)
    push!(indexer.fIy, var.factor)

end

function populate!(indexer::FrameOptIndexer, var::CoupledVariable{SectionVariable{SectionJ}})

    push!(indexer.iJ, var.i)
    push!(indexer.iJg, var.iglobal)
    push!(indexer.fJ, var.factor)
    
end