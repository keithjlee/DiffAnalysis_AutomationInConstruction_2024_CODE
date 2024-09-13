"""
    TrussResults

Results of a differentiable structural analysis of a truss model.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- A: Area of all elements
- L: Length of all elements
- K: Elemental stiffness matrices in GCS
- R: Elemental transformation matrices
- U: Displacement vector of all nodes
"""
struct TrussResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    A::Vector{Float64}
    L::Vector{Float64}
    K::Vector{Matrix{Float64}}
    R::Vector{Matrix{Float64}}
    U::Vector{Float64}
end

"""
    NetworkResults

Results of a differentiable FDM analysis of a network.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- Q: force densities of all elements
"""
struct NetworkResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    Q::Vector{Float64}
end

"""
    TrussResults

Results of a differentiable structural analysis of a truss model.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- A: Area of all elements
- L: Length of all elements
- K: Elemental stiffness matrices in GCS
- R: Elemental transformation matrices
- U: Displacement vector of all nodes
"""
struct FrameResults
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    A::Vector{Float64}
    Ix::Vector{Float64}
    Iy::Vector{Float64}
    J::Vector{Float64}
    L::Vector{Float64}
    K::Vector{Matrix{Float64}}
    R::Vector{Matrix{Float64}}
    U::Vector{Float64}
end

"""
    GeometricProperties

Get the primary and secondary geometric properties of a truss structure in a differentiable format.

Fields:
- X: x-position of all nodes
- Y: y-position of all nodes
- Z: z-position of all nodes
- evecs: [nel × 3] matrix where row_i = nodes[iend].position - nodes[istart].position
- A: Area of all elements
- L: Length of all elements
- 
"""
struct GeometricProperties
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
    evecs::Matrix{Float64}
    A::Vector{Float64}
    L::Vector{Float64}

    GeometricProperties(x, y, z, e, a, l) = new(x, y, z, e, a, l)

    function GeometricProperties(design_variables::Vector{Float64}, opt_params::TrussOptParams)

        #populate values
        X = opt_params.indexer.activeX ? add_values(opt_params.X, opt_params.indexer.iX, design_variables[opt_params.indexer.iXg] .* opt_params.indexer.fX) : opt_params.X
        Y = opt_params.indexer.activeY ? add_values(opt_params.Y, opt_params.indexer.iY, design_variables[opt_params.indexer.iYg] .* opt_params.indexer.fY) : opt_params.Y
        Z = opt_params.indexer.activeZ ? add_values(opt_params.Z, opt_params.indexer.iZ, design_variables[opt_params.indexer.iZg] .* opt_params.indexer.fZ) : opt_params.Z
        A = opt_params.indexer.activeA ? replace_values(opt_params.A, opt_params.indexer.iA, design_variables[opt_params.indexer.iAg] .* opt_params.indexer.fA) : opt_params.A

        # [nₑₗ × 3] matrix where row i is the vector representation of element i, from the start node to the end node; ||vecₑ|| = Lₑ
        evecs = get_element_vectors(X, Y, Z, opt_params)
        L = get_element_lengths(evecs)

        return new(
            X,
            Y,
            Z,
            evecs,
            A,
            L
        )

    end

    function GeometricProperties(design_variables::Vector{Float64}, opt_params::FrameOptParams)

        #populate values
        X = opt_params.indexer.activeX ? add_values(opt_params.X, opt_params.indexer.iX, design_variables[opt_params.indexer.iXg] .* opt_params.indexer.fX) : opt_params.X
        Y = opt_params.indexer.activeY ? add_values(opt_params.Y, opt_params.indexer.iY, design_variables[opt_params.indexer.iYg] .* opt_params.indexer.fY) : opt_params.Y
        Z = opt_params.indexer.activeZ ? add_values(opt_params.Z, opt_params.indexer.iZ, design_variables[opt_params.indexer.iZg] .* opt_params.indexer.fZ) : opt_params.Z
        A = opt_params.indexer.activeA ? replace_values(opt_params.A, opt_params.indexer.iA, design_variables[opt_params.indexer.iAg] .* opt_params.indexer.fA) : opt_params.A

        # [nₑₗ × 3] matrix where row i is the vector representation of element i, from the start node to the end node; ||vecₑ|| = Lₑ
        evecs = get_element_vectors(X, Y, Z, opt_params)
        L = get_element_lengths(evecs)

        return new(
            X,
            Y,
            Z,
            evecs,
            A,
            L
        )

    end
end

export GeometricProperties_noadjoint
function GeometricProperties_noadjoint(design_variables::Vector{Float64}, opt_params::TrussOptParams)
    #populate values
    X = opt_params.indexer.activeX ? add_values(opt_params.X, opt_params.indexer.iX, design_variables[opt_params.indexer.iXg] .* opt_params.indexer.fX) : opt_params.X
    Y = opt_params.indexer.activeY ? add_values(opt_params.Y, opt_params.indexer.iY, design_variables[opt_params.indexer.iYg] .* opt_params.indexer.fY) : opt_params.Y
    Z = opt_params.indexer.activeZ ? add_values(opt_params.Z, opt_params.indexer.iZ, design_variables[opt_params.indexer.iZg] .* opt_params.indexer.fZ) : opt_params.Z
    A = opt_params.indexer.activeA ? replace_values(opt_params.A, opt_params.indexer.iA, design_variables[opt_params.indexer.iAg] .* opt_params.indexer.fA) : opt_params.A

    # [nₑₗ × 3] matrix where row i is the vector representation of element i, from the start node to the end node; ||vecₑ|| = Lₑ
    evecs = get_element_vectors_noadjoint(X, Y, Z, opt_params)
    L = get_element_lengths_noadjoint(evecs)

    return GeometricProperties(
            X,
            Y,
            Z,
            evecs,
            A,
            L
        )
end