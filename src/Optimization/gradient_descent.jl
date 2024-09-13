struct GDresults
    minimum::Float64
    minimizer::Vector{Float64}
    n_iter::Int64
    stop_criteria::Symbol
    loss_history::Vector{Float64}
    x_history::Vector{Vector{Float64}}
    ftol_rel::Float64
    ftol_abs::Float64
    step_size::Float64
    max_iter::Int64
    normalize_grad::Bool
end

function gradient_descent(fn::Function, x0::Vector{Float64}; ftol_rel = 1e-4, ftol_abs = 1e-6, step_size = 0.1, max_iter = 300, normalize_grad = true)

    #initialize
    loss_store = Vector{Float64}()
    x_store = Vector{Vector{Float64}}()
    stop_criteria = :init
    
    iter = 1
    reldiff = Inf
    adiff = Inf

    x = deepcopy(x0)
    while iter ≤ max_iter
        val, g = withgradient(fn, x)

        #define gradient
        grad = normalize_grad ? g[1] ./ maximum(abs.(g[1])) : g[1]

        #save history
        push!(loss_store, val)
        push!(x_store, x)


        #check
        if iter > 1
            adiff = abs(val - loss_store[end-1])
            reldiff = adiff / loss_store[end-1]
        end

        if adiff < ftol_abs
            stop_criteria = :ABSTOL_REACHED
            break
        elseif reldiff < ftol_rel
            stop_criteria = :RELTOL_REACHED
            break
        end

        #update
        x -= grad .* step_size
        iter += 1
    end

    # update stop criteria if max iters reached
    stop_criteria == :init && (stop_criteria = :MAXITER_REACHED)

    return GDresults(last(loss_store), last(x_store), iter, stop_criteria, loss_store, x_store, ftol_rel, ftol_abs, step_size, max_iter, normalize_grad)

end