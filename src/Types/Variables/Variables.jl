include("SpatialVariables.jl")
include("AreaVariables.jl")
include("NetworkVariables.jl")
include("FrameVariables.jl")
include("CoupledVariables.jl")

const TrussVariable = Union{SpatialVariable, AreaVariable, CoupledVariable, VectorSpatialVariable}
const NetworkVariable = Union{SpatialVariable, QVariable, CoupledVariable, VectorSpatialVariable}
const FrameVariable = Union{SpatialVariable, AreaVariable, SectionVariable, CoupledVariable, VectorSpatialVariable}

function process_variables!(vars::Vector{T}) where {T<:AbstractVariable}

    # get the variable type indices
    i_independent = findall(typeof.(vars) .<: IndependentVariable)
    i_coupled = findall(typeof.(vars) .<: CoupledVariable)

    # number of independent variables
    n_independent = length(i_independent)

    # collectors
    uid_to_globalid = Dict{UInt64, Int64}()
    vals = Vector{Float64}(undef, n_independent)
    lb = Vector{Float64}(undef, n_independent)
    ub = Vector{Float64}(undef, n_independent)

    # populate independent variable data
    i_global = 1
    for (i, var) in enumerate(vars[i_independent])
    
        var.iglobal = i_global
        uid_to_globalid[objectid(var)] = i_global
        vals[i] = var.val
        lb[i] = var.lb
        ub[i] = var.ub
        i_global += 1
        
    end

    #populate coupled variables
    if length(i_coupled) > 0
        for var in vars[i_coupled]
            var.iglobal = uid_to_globalid[var.target]
        end
    end

    return vals, lb, ub
end