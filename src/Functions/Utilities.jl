"""
    replace_values(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

Replace the values of `values[indices]` with the values in `newvalues` in a differentiable way. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function replace_values(values::Vector{Float64}, indices::Vector{Int64}, newvalues)
    
    v2 = zero(values) + values
    v2[indices] .= newvalues

    return v2

end

"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement.

df/dnewvalues = df/dvalues ⋅ dvalues / dnewvalues = v̄ ⋅ dvalues / dnewvalues

Where v̄ = [nvalues × 1], dvalues/dnewvalues = [nvalues × nnewvalues] so:

df/dnewvalues = (dvalues/dnewvalues)ᵀv̄ = [nnewvalues × 1]

Is simply the values of v̄ at the indices of the new values.
"""
function ChainRulesCore.rrule(::typeof(replace_values), values::Vector{Float64}, indices::Vector{Int64}, newvalues)

    v = replace_values(values, indices, newvalues)

    function replace_values_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, replace_values_pullback 
end

"""
    replace_values_buffer(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

Return a new vector `x` where `x[indices] = newvalues` and `x[.!indices] == values[.!indices]` using Zygote's buffer system such that this function is differentiable with respect to values.
"""
function replace_values_buffer(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

    buf = Zygote.bufferfrom(values)

    buf[indices] = newvalues

    return copy(buf)

end

"""
    add_values(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

Add the values of `increments` to the current values in `values` at `indices`. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function add_values(values::Vector{Float64}, indices::Vector{Int64}, increments)

    v2 = zero(values) + values
    @views v2[indices] .+= increments

    return v2
end


"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement
"""
function ChainRulesCore.rrule(::typeof(add_values), values::Vector{Float64}, indices::Vector{Int64}, increments)

    v = add_values(values, indices, increments)

    function add_values_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, add_values_pullback 
end

"""
    add_values_buffer(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

Return a new vector `x` where `x[indices] = newvalues + values[indices]` and `x[.!indices] == values[.!indices]` using Zygote's buffer system such that this function is differentiable with respect to values.
"""
function add_values_buffer(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    buf = Zygote.Buffer(values)
    copyto!(buf, values)

    buf[indices] += increments

    return copy(buf)
end