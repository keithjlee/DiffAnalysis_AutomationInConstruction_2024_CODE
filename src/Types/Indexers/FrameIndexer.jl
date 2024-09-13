mutable struct FrameOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    fX::Vector{Float64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    fY::Vector{Float64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    fZ::Vector{Float64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
    fA::Vector{Float64}
    iIx::Vector{Int64}
    iIxg::Vector{Int64}
    fIx::Vector{Float64}
    iIy::Vector{Int64}
    iIyg::Vector{Int64}
    fIy::Vector{Float64}
    iJ::Vector{Int64}
    iJg::Vector{Int64}
    fJ::Vector{Float64}
    activeX::Bool
    activeY::Bool
    activeZ::Bool
    activeA::Bool
    activeIx::Bool
    activeIy::Bool
    activeJ::Bool
end

function populate!(indexer::FrameOptIndexer, var::SectionVariable{SectionA})
    push!(indexer.iA, var.i)
    push!(indexer.iAg, var.iglobal)
    push!(indexer.fA, 1.0)

    indexer.activeA = true
end

function populate!(indexer::FrameOptIndexer, var::SectionVariable{SectionIx})
    push!(indexer.iIx, var.i)
    push!(indexer.iIxg, var.iglobal)
    push!(indexer.fIx, 1.0)

    indexer.activeIx = true
end

function populate!(indexer::FrameOptIndexer, var::SectionVariable{SectionIy})
    push!(indexer.iIy, var.i)
    push!(indexer.iIyg, var.iglobal)
    push!(indexer.fIy, 1.0)

    indexer.activeIy = true
end

function populate!(indexer::FrameOptIndexer, var::SectionVariable{SectionJ})
    push!(indexer.iJ, var.i)
    push!(indexer.iJg, var.iglobal)
    push!(indexer.fJ, 1.0)

    indexer.activeJ = true
end


function FrameOptIndexer(vars::Vector{T}) where T<:FrameVariable
    indexer = FrameOptIndexer(
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        false,
        false,
        false,
        false,
        false,
        false,
        false
        )

    # iX = Vector{Int64}()
    # iXg = Vector{Int64}()
    # fX = Vector{Float64}()
    # iY = Vector{Int64}()
    # iYg = Vector{Int64}()
    # fY = Vector{Float64}()
    # iZ = Vector{Int64}()
    # iZg = Vector{Int64}()
    # fZ = Vector{Float64}()
    # iA = Vector{Int64}()
    # iAg = Vector{Int64}()
    # fA = Vector{Float64}()
    # iIx = Vector{Int64}()
    # iIxg = Vector{Int64}()
    # fIx = Vector{Float64}()
    # iIy = Vector{Int64}()
    # iIyg = Vector{Int64}()
    # fIy = Vector{Float64}()
    # iJ = Vector{Int64}()
    # iJg = Vector{Int64}()
    # fJ = Vector{Float64}()
    # activeX = false
    # activeY = false
    # activeZ = false
    # activeA = false
    # activeIx = false
    # activeIy = false
    # activeJ = false

    for var in vars
        populate!(indexer, var)
    end

    indexer
end