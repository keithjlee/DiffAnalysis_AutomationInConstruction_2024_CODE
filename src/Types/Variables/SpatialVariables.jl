abstract type SpatialVarDirections end
struct SpatialX <: SpatialVarDirections end
struct SpatialY <: SpatialVarDirections end
struct SpatialZ <: SpatialVarDirections end

"""
    SpatialVariable <: IndependentVariable

A variable tied to the spatial position of a node in a single axis.

```julia
SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
SpatialVariable(node::Union{Asap.AbstractNode, Asap.FDMnode}, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
```
"""
mutable struct SpatialVariable{T<:SpatialVarDirections} <: IndependentVariable
    i::Int64 #index of node, e.g. X[i] is the spatial variable
    val::Float64 #value
    lb::Float64 #lower bound of variable
    ub::Float64 #upper bound of variable
    iglobal::Int64 # position in the vector of active design variables
end

axis_to_spatial_type = Dict(
    :x => SpatialX,
    :X => SpatialX,
    :y => SpatialY,
    :Y => SpatialY,
    :z => SpatialZ,
    :Z => SpatialZ
)

function SpatialVariable(node::Asap.AbstractNode, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
    #convert axis to relevant type
    T = axis_to_spatial_type[axis]
    return SpatialVariable{T}(node.nodeID, value, lowerbound, upperbound, 0)
end

function SpatialVariable(index::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)
    #convert axis to relevant type
    T = axis_to_spatial_type[axis]
    return SpatialVariable{T}(index, value, lowerbound, upperbound, 0)
end

function SpatialVariable(node::Asap.FDMnode, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

    @assert node.dof == false "FDM spatial variables only apply to anchor (fixed) nodes"

    #convert axis to relevant type
    T = axis_to_spatial_type[axis]
    return SpatialVariable{T}(node.nodeID, value, lowerbound, upperbound, 0)
end

SpatialVariable(node::Union{Asap.Node, Asap.TrussNode, Asap.FDMnode}, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z) = SpatialVariable(node, 0.0, lowerbound, upperbound, axis)

mutable struct VectorSpatialVariable <: IndependentVariable
    i::Int64 #index of node
    val::Float64 #value
    vec::Vector{Float64} #directional vector of variable. gives pos = pos₀ + vec * val
    lb::Float64 #lower bound of val
    ub::Float64 #upper bound of val
    iglobal::Int64 #position in the vector of active design variables

    function VectorSpatialVariable(node::Asap.AbstractNode, vector::Vector{Float64}, value::Float64, lowerbound::Float64, upperbound::Float64)

        @assert length(vector) == 3

        return new(node.nodeID, value, vector, lowerbound, upperbound, 0)

    end

end
