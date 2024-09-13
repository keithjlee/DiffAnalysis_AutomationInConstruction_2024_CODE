function disp_stress_cstr(U::Vector{Float64}, Kes::Vector{Matrix{Float64}}, Aes::Vector{Float64}, p::P, dstart::Int64, dinc::Int64, dmax::Float64, fmax::Float64) where {P<:AbstractOptParams}

    disps = abs.(U[dstart:dinc:end])

    stresses = Vector{Float64}(undef, length(Aes))
    for i in eachindex(stresses)
        @views stresses[i] = norm(Kes[i][4:6,:] * U[p.dofids[i]]) / Aes[i]
    end

    return [
        (disps .- dmax);
        (stresses .- fmax);
    ]

end

function ChainRulesCore.rrule(::typeof(disp_stress_cstr), U::Vector{Float64}, Kes::Vector{Matrix{Float64}}, Aes::Vector{Float64}, p::P, dstart::Int64, dinc::Int64, dmax::Float64, fmax::Float64) where {P<:AbstractOptParams}


    disps = abs.(U[dstart:dinc:end])
    forces = [norm(k[4:6,:] * U[id]) for (k, id) in zip(Kes, p.dofids)]

    out = [
        (disps .- dmax);
        (forces ./ Aes .- fmax);
    ]

    disp_indices = 1:length(disps)
    stress_indices = length(disps)+1:length(out)

    function _pullback(C̄)

        dU = zero(U)
        dKs = zero.(Kes)
        dAs = -forces ./ Aes.^2 .* C̄[stress_indices]

        #displacement pullbacks
        @views dU[dstart:dinc:end] .+= sign.(U[dstart:dinc:end]) .* C̄[disp_indices]

        #stress pullbacks
        for i in eachindex(dKs)
            cbar = C̄[stress_indices[i]]
            Ksub = Kes[i][4:6, :]
            uk = U[p.dofids[i]]

            # if forces[i] > 1e-6
                dKs[i][4:6, :] .= 1 / forces[i] / Aes[i] * Ksub * uk * uk' .* cbar
                dU[p.dofids[i]] .+= cbar / forces[i] / Aes[i] * Ksub' * Ksub * uk
            # end

        end

        return (NoTangent(), dU, dKs, dAs, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())

    end

    return out, _pullback

end

function stress_cstr(U::Vector{Float64}, Kes::Vector{Matrix{Float64}}, Aes::Vector{Float64}, p::TrussOptParams, fmax::Float64)

    stresses = Vector{Float64}(undef, length(Aes))
    for i in eachindex(stresses)
        @views stresses[i] = norm(Kes[i][4:6,:] * U[p.dofids[i]]) / Aes[i]
    end

    return stresses .- fmax
end

function ChainRulesCore.rrule(::typeof(stress_cstr), U::Vector{Float64}, Kes::Vector{Matrix{Float64}}, Aes::Vector{Float64}, p::TrussOptParams, fmax::Float64)

    forces = Vector{Float64}(undef, length(Aes))
    for i in eachindex(forces)
        @views forces[i] = norm(Kes[i][4:6,:] * U[p.dofids[i]])
    end

    function stress_cstr_pullback(C̄)

        dAs = -forces ./ Aes.^2 .* C̄

    end

end