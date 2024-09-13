"""
    get_element_vectors(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::TrussOptParams)

Get the [nₑ × 3] matrix where each row is the [x,y,z] vector of an element
"""
function get_element_vectors(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::AbstractOptParams)
    p.C * [X Y Z]
end

function get_element_vectors_noadjoint(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::AbstractOptParams)
    p.C * [X Y Z]
end

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
function ChainRulesCore.rrule(::typeof(get_element_vectors), X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, p::AbstractOptParams)
    v = get_element_vectors(X, Y, Z, p)

    function get_element_vectors_pullback(v̄)
        
        dv = p.C' * v̄

        return (NoTangent(), dv[:,1], dv[:,2], dv[:,3], NoTangent())
        
    end

    return v, get_element_vectors_pullback
end


"""
    get_element_lengths(XYZ::Matrix{Float64})

Get the [nₑ × 1] vector of element lengths
"""
function get_element_lengths(XYZ::Matrix{Float64})
    norm.(eachrow(XYZ))
end

function get_element_lengths_noadjoint(XYZ::Matrix{Float64})
    norm.(eachrow(XYZ))
end

"""
For a single element, given its vector representation:

L(element) = ||xyzₑ|| = √(x² + y² + z²)
dL/dx = x/L

dg/dx = dg/dL ⋅ dL/dx = L̄ dL/dx
"""
function ChainRulesCore.rrule(::typeof(get_element_lengths), XYZ::Matrix{Float64})
    l = get_element_lengths(XYZ)

    function l_pullback(l̄)
        dl = l̄ ./ l .* XYZ 
        
        return (NoTangent(), dl)
    end

    return l, l_pullback
end

"""
    get_normalized_element_vectors(XYZ::Matrix{Float64}, Ls::Vector{Float64})

Get the unit vector representation of elements (local x axis)
"""
function get_normalized_element_vectors(XYZ::Matrix{Float64}, Ls::Vector{Float64})
    XYZ ./ Ls
end

function get_normalized_element_vectors_noadjoint(XYZ::Matrix{Float64}, Ls::Vector{Float64})
    XYZ ./ Ls
end

"""
g = f(XYZn)

dg/dXYZ = df/dXYZn ⋅ dXYZn/dXYZ = v̄ ⋅ dXYZn/dXYZ = v̄ ⋅ [1 1 1; 1 1 1; ...] ./ L
dg/dL = v̄ ⋅ dXYZn/dL = -v̄ ⋅ XYZ / L^2 = -v̄ ⋅ XYZn / L
"""
function ChainRulesCore.rrule(::typeof(get_normalized_element_vectors), XYZ::Matrix{Float64}, Ls::Vector{Float64})
    XYZn = get_normalized_element_vectors(XYZ, Ls)

    function getnormv_pullback(v̄)
        dxyz = v̄ ./ Ls .* ones(size(XYZn)...)

        dL = sum.(eachrow(-XYZn .* v̄ ./ Ls))

        return (NoTangent(), dxyz, dL)
    end

    return XYZn, getnormv_pullback
end