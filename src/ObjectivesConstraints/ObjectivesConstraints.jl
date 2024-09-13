function objective_volume(x::Vector{Float64}, p::TrussOptParams)

    geo = GeometricProperties(x, p)

    return dot(geo.L, geo.A)

end

function objective_compliance(x::Vector{Float64}, p::TrussOptParams; alg = UMFPACKFactorization())

    res = solve_truss(x, p; linsolve_alg = alg)

    return dot(res.U, p.P)
end

function objective_compliance(x::Vector{Float64}, p::FrameOptParams; alg = UMFPACKFactorization())

    res = solve_frame(x, p; linsolve_alg = alg)

    return dot(res.U, p.P)
end

function constraint_volume(x::Vector{Float64}, p::Union{TrussOptParams, FrameOptParams}, Vmax::Float64)

    geo = GeometricProperties(x, p)

    return dot(geo.A, geo.L) - Vmax

end

function constraint_length(x::Vector{Float64}, p::Union{TrussOptParams, FrameOptParams}, Lmax::Float64)

    geo = GeometricProperties(x, p)

    return geo.L .- Lmax

end