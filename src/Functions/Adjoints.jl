"""
The elemental vectors are derived from V = C * XYZ where:

- C: [nₑ × nₙ] matrix defining the topology of the structure
- XYZ: [nₙ × 3] matrix defining the positions of nodes

Given an downstream function g = f(V), then the gradient of g w/r/t an input argument (e.g., X) is:

dg/dX = df/dV ⋅ dV/dX = V̄ ⋅ dV/dX

And dV/dX = d/dX (C X) = C

Such that:

dg/dX = CᵀV̄

and likewise for Y, Z.
"""
function ChainRulesCore.rrule(::typeof(getevecs), X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::AbstractOptParams)
    v = getevecs(X, Y, Z, p)

    function getevecs_pullback(v̄)
        
        dv = p.C' * v̄

        return (NoTangent(), dv[:,1], dv[:,2], dv[:,3], NoTangent())
        
    end

    return v, getevecs_pullback
end

"""
For a single element, given its vector representation:

L(element) = ||xyzₑ|| = √(x² + y² + z²)
dL/dx = x/L

dg/dx = dg/dL ⋅ dL/dx = L̄ dL/dx
"""
function ChainRulesCore.rrule(::typeof(getlengths), XYZ::Matrix{Float64})
    l = getlengths(XYZ)

    function l_pullback(l̄)
        dl = l̄ ./ l .* XYZ 
        
        return (NoTangent(), dl)
    end

    return l, l_pullback
end

"""
g = f(XYZn)

dg/dXYZ = df/dXYZn ⋅ dXYZn/dXYZ = v̄ ⋅ dXYZn/dXYZ = v̄ ⋅ [1 1 1; 1 1 1; ...] ./ L
dg/dL = v̄ ⋅ dXYZn/dL = -v̄ ⋅ XYZ / L^2 = -v̄ ⋅ XYZn / L
"""
function ChainRulesCore.rrule(::typeof(getnormalizedevecs), XYZ::Matrix{Float64}, Ls::Vector{Float64})
    XYZn = getnormalizedevecs(XYZ, Ls)

    function getnormv_pullback(v̄)
        dxyz = v̄ ./ Ls .* ones(size(XYZn)...)

        dL = sum.(eachrow(-XYZn .* v̄ ./ Ls))

        return (NoTangent(), dxyz, dL)
    end

    return XYZn, getnormv_pullback
end

"""
For a single element:
dKe/dR = 2KeΓ
dKe/dk = ΓkΓᵀ
"""
function ChainRulesCore.rrule(::typeof(getglobalks), rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})
        
    kgs = getglobalks(rs, ks)

    function kgs_pullback(k̄)
        dr = (2 .* ks .* rs) .* k̄
        dk = rs .* k̄ .* transpose.(rs)

        return (NoTangent(), dr, dk)
    end

    return kgs, kgs_pullback
end

"""
Adjoint w/r/t element variables E, A, L for the local stiffness matrix of a truss element
"""
function ChainRulesCore.rrule(::typeof(ktruss), E::Float64, A::Float64, L::Float64)
    k = ktruss(E, A, L)

    function ktruss_pullback(k̄)

        # ∇E = dot(k̄, (A / L * [1 -1; -1 1]))
        ∇A = dot(k̄, (E / L * [1 -1; -1 1]))
        ∇L = dot(k̄, (- E * A / L^2 * [1 -1; -1 1]))

        return (NoTangent(), NoTangent(), ∇A, ∇L)
    end

    return k, ktruss_pullback
end

dRdx = [1. 0. 0. 0. 0. 0.; 0. 0. 0. 1. 0. 0.]
dRdy = [0. 1. 0. 0. 0. 0.; 0. 0. 0. 0. 1. 0.]
dRdz = [0. 0. 1. 0. 0. 0.; 0. 0. 0. 0. 0. 1.]

"""
Adjoint for global transformation matrix
"""
function ChainRulesCore.rrule(::typeof(Rtruss), Cx::Float64, Cy::Float64, Cz::Float64)
    R = Rtruss(Cx, Cy, Cz)

    function Rtruss_pullback(R̄)
        return (NoTangent(),
            dot(R̄, dRdx),
            dot(R̄, dRdy),
            dot(R̄, dRdz))
    end

    return R, Rtruss_pullback
end

function ChainRulesCore.rrule(::typeof(Rtruss), Cxyz::SubArray)
    R = Rtruss(Cxyz)

    function Rtruss_pullback(R̄)
        return (NoTangent(),
            [dot(R̄, dRdx),
            dot(R̄, dRdy),
            dot(R̄, dRdz)])
    end

    return R, Rtruss_pullback
end

function ChainRulesCore.rrule(::typeof(getRmatrices), XYZn::Matrix{Float64})
    rs = getRmatrices(XYZn)

    function getRmatrices_pullback(R̄)
        dRdC = zero(XYZn)

        for i in axes(R̄, 1)
            dRdC[i, 1] = dot(R̄[i], dRdx)
            dRdC[i, 2] = dot(R̄[i], dRdy)
            dRdC[i, 3] = dot(R̄[i], dRdz)
        end

        return (NoTangent(), dRdC)
    end

    return rs, getRmatrices_pullback
end

"""
The sensitivity of K w/r/t an elemental Kₑ is the proportional stiffness added to K from Kₑ

Output is a vector: [nElements × [nDOFe × nDOFe]] of elemental stiffness matrix sensitivites
"""
function ChainRulesCore.rrule(::typeof(assembleglobalK), Eks::Vector{Matrix{Float64}}, p::TrussOptParams)
    K = assembleglobalK(Eks, p)

    function K_pullback(K̄)
        dEks = [Matrix(K̄[id, id]) for id in p.dofids]

        return NoTangent(), dEks, NoTangent()
    end

    return K, K_pullback
end

"""
u = inv(K) * P

if obj = f(u), then the gradient of obj with respect to an independent variable x is achieved through the chain rule:

dObj/dx = df/du ⋅ du/dK ⋅ ... = ū ⋅ dy/dK ⋅ ...

For this rule, we are concerned with finding du/dK, or the [ndof × ndof] matrix of sensitivites that we can propagate backwards to the final objective.

Given df/du = [ndof × 1] = ū is the gradient of the objective function with respect to displacements u, the sensitivity is:

du/dK = - uᵀ ⊗ K⁻¹
df/dK = du/dK ū = - (uᵀ ⊗ K⁻¹)ū

Can be rearranged such that:
ΔK = K⁻¹ū

df/dK = -uᵀ ⊗ ΔK

Which is an [ndof × ndof] matrix where:

Columnᵢ = uᵢ .* ΔK
"""
function ChainRulesCore.rrule(::typeof(solveU), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    u = solveU(K, p)

    function solveU_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = cg(K[p.freeids, p.freeids], ū)

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solveU_pullback
end

"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement.

df/dnewvalues = df/dvalues ⋅ dvalues / dnewvalues = v̄ ⋅ dvalues / dnewvalues

Where v̄ = [nvalues × 1], dvalues/dnewvalues = [nvalues × nnewvalues] so:

df/dnewvalues = (dvalues/dnewvalues)ᵀv̄ = [nnewvalues × 1]

Is simply the values of v̄ at the indices of the new values.
"""
function ChainRulesCore.rrule(::typeof(replacevalues), values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

    v = replacevalues(values, indices, newvalues)

    function replacevalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, replacevalues_pullback 
end


"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement
"""
function ChainRulesCore.rrule(::typeof(addvalues), values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    v = addvalues(values, indices, increments)

    function addvalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, addvalues_pullback 
end