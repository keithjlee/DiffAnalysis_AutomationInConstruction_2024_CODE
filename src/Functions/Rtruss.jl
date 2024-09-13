const dRdx_truss = [1. 0. 0. 0. 0. 0.; 0. 0. 0. 1. 0. 0.]
const dRdy_truss = [0. 1. 0. 0. 0. 0.; 0. 0. 0. 0. 1. 0.]
const dRdz_truss = [0. 0. 1. 0. 0. 0.; 0. 0. 0. 0. 0. 1.]

"""
    r_truss(Cx::Float64, Cy::Float64, Cz::Float64)

Transformation matrix for truss element
"""
function r_truss(Cx::Float64, Cy::Float64, Cz::Float64)
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

function r_truss_noadjoint(Cx::Float64, Cy::Float64, Cz::Float64)
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

"""
Adjoint for global transformation matrix
"""
function ChainRulesCore.rrule(::typeof(r_truss), Cx::Float64, Cy::Float64, Cz::Float64)
    R = r_truss(Cx, Cy, Cz)

    function r_truss_pullback(R̄)
        return (NoTangent(),
            dot(R̄, dRdx_truss),
            dot(R̄, dRdy_truss),
            dot(R̄, dRdz_truss))
    end

    return R, r_truss_pullback
end

"""
    r_truss(Cxyz::SubArray)

Transformation of sliced view of matrix of element vectors
"""
function r_truss(Cxyz::SubArray)
    Cx, Cy, Cz = Cxyz
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

function r_truss_noadjoint(Cxyz::SubArray)
    Cx, Cy, Cz = Cxyz
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

function ChainRulesCore.rrule(::typeof(r_truss), Cxyz::SubArray)
    R = r_truss(Cxyz)

    function r_truss_pullback(R̄)
        return (NoTangent(),
            [dot(R̄, dRdx_truss),
            dot(R̄, dRdy_truss),
            dot(R̄, dRdz_truss)])
    end

    return R, r_truss_pullback
end

"""
    r_truss(XYZn::Matrix{Float64})

Get all transformation matrices from a [nₑ × 3] matrix of all normalized element local x vectors
"""
function r_truss(XYZn::Matrix{Float64})
    r_truss.(eachrow(XYZn))
end

function r_truss_noadjoint(XYZn::Matrix{Float64})
    r_truss_noadjoint.(eachrow(XYZn))
end

function ChainRulesCore.rrule(::typeof(r_truss), XYZn::Matrix{Float64})
    rs = r_truss(XYZn)

    function r_truss_pullback(R̄)
        dRdC = zero(XYZn)

        for i in axes(R̄, 1)
            dRdC[i, 1] = dot(R̄[i], dRdx_truss)
            dRdC[i, 2] = dot(R̄[i], dRdy_truss)
            dRdC[i, 3] = dot(R̄[i], dRdz_truss)
        end

        return (NoTangent(), dRdC)
    end

    return rs, r_truss_pullback
end