"""
    CoupledVariable <: AbstractVariable

A variable whose value points to an existing variable.

```julia
var = CoupledVariable(variable::Union{TrussNode, TrussElement, FDMelement}, reference::Union{SpatialVariable, AreaVariable, QVaraible}, factor = 1.0)
```

Where `factor` is a scalar factor applied to the value of `reference` before assigning to `variable`. IE to enforce a mirrored value, use `factor = -1.0`.
"""
mutable struct CoupledVariable{T<:IndependentVariable} <: AbstractVariable
    i::Int64
    target::UInt64
    factor::Union{Float64, Vector{Float64}}
    iglobal::Int64
end

function CoupledVariable(node::Asap.AbstractNode, ref::T, factor::Float64 = 1.0) where {T<:SpatialVariable}

    CoupledVariable{T}(node.nodeID, objectid(ref), factor, 0)
end

function CoupledVariable(index::Int64, ref::T, factor::Float64 = 1.0) where {T<:SpatialVariable}

    CoupledVariable{T}(index, objectid(ref), factor, 0)
end


function CoupledVariable(element::Asap.AbstractElement, ref::AreaVariable, factor::Float64 = 1.0)

    @assert factor > 0

    CoupledVariable{AreaVariable}(element.elementID, objectid(ref), factor, 0)
end

function CoupledVariable(element::Asap.Element, ref::T, factor::Float64 = 1.0) where {T<:SectionVariable}

    @assert factor > 0

    CoupledVariable{T}(element.elementID, objectid(ref), factor, 0)
end

function CoupledVariable(element::Asap.FDMelement, ref::QVariable, factor::Float64 = 1.0)

    @assert factor > 0

    CoupledVariable{QVariable}(element.elementID, objectid(ref), factor, 0)
end

function CoupledVariable(node::Asap.AbstractNode, ref::VectorSpatialVariable, vec::Vector{Float64})

    @assert length(vec) == 3

    CoupledVariable{VectorSpatialVariable}(node.nodeID, objectid(ref), vec, 0)

end

function CoupledVariable(node::Asap.AbstractNode, ref::VectorSpatialVariable, factor::Union{Float64, Vector{Float64}} = 1.0)

    @assert length(factor) == 3 || length(factor) == 1

    CoupledVariable{VectorSpatialVariable}(node.nodeID, objectid(ref), ref.vec .* factor, 0)

end