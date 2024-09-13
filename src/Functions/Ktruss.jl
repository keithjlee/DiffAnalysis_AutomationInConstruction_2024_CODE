"""
    k_truss(E::Float64, A::Float64, L::Float64)

Get the element truss stiffness matrix in the local coordinate system
"""
function k_truss(E::Float64, A::Float64, L::Float64)
    E * A / L * [1 -1; -1 1]
end

function k_truss_noadjoint(E::Float64, A::Float64, L::Float64)
    E * A / L * [1 -1; -1 1]
end

"""
Adjoint w/r/t element variables E, A, L for the local stiffness matrix of a truss element
"""
function ChainRulesCore.rrule(::typeof(k_truss), E::Float64, A::Float64, L::Float64)
    k = k_truss(E, A, L)

    function k_truss_pullback(k̄)

        ∇E = dot(k̄, (A / L * [1 -1; -1 1]))
        ∇A = dot(k̄, (E / L * [1 -1; -1 1]))
        ∇L = dot(k̄, (- E * A / L^2 * [1 -1; -1 1]))

        return (NoTangent(), ∇E, ∇A, ∇L)
    end

    return k, k_truss_pullback
end

"""
    getlocalks(E::Vector{Float64}, A::Vector{Float64}, L::Vector{Float64})

Broadcast function for vectors of element properties
"""
function getlocalks(E::Vector{Float64}, A::Vector{Float64}, L::Vector{Float64})
    k_truss.(E, A, L)
end

function getlocalks_noadjoint(E::Vector{Float64}, A::Vector{Float64}, L::Vector{Float64})
    k_truss_noadjoint.(E, A, L)
end